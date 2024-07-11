!> QUINCY calculate vegetation dynamics
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
!>#### calculate vegetation dynamics, e.g., establishment and mortality
!>
MODULE mo_q_veg_dynamics
#ifndef __NO_QUINCY__

  USE mo_kind,                 ONLY: wp
  USE mo_jsb_impl_constants,   ONLY: true, false, test_false_true
  USE mo_jsb_control,          ONLY: debug_on
  USE mo_exception,            ONLY: message
  USE mo_jsb_math_constants,   ONLY: one_year, one_day, dtime, eps8, eps4

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID, VEG_BGCM_LITTERFALL_ID, VEG_BGCM_ESTABLISHMENT_ID

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: update_veg_dynamics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_veg_dynamics'

CONTAINS

  ! ======================================================================================================= !
  !>Calculate establishment and mortality of a population
  !>
  !>  (currently static mortality rate only)
  !>  a background_mort_rate exists for each trees and grasses, defined as mortality rate per year
  !>
  SUBROUTINE update_veg_dynamics(tile, options)

    USE mo_jsb_lctlib_class,      ONLY: t_lctlib_element
    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: VEG_, Q_RAD_, Q_PHENO_, SPQ_
    USE mo_jsb4_forcing,          ONLY: forcing_options_quincy
    USE mo_veg_constants,         ONLY: ITREE
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_config(VEG_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_RAD_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),        POINTER :: model                  !< the model
    TYPE(t_lnd_bgcm_store),   POINTER :: bgcm_store             !< the bgcm store of this tile
    TYPE(t_lctlib_element),   POINTER :: lctlib                 !< land-cover-type library - parameter across pft's
    REAL(wp), DIMENSION(options%nc)   :: fpc                    !< current foliage projective cover
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_veg_dynamics'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    dsl4jsb_Def_mt2L2D :: veg_litterfall_mt
    dsl4jsb_Def_mt1L2D :: veg_establishment_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(Q_RAD_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_RAD_ 2D
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! SPQ_ 2D
    dsl4jsb_Real2D_onChunk      :: w_soil_root_theta
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: do_cohort_harvest
    dsl4jsb_Real2D_onChunk      :: dens_ind
    dsl4jsb_Real2D_onChunk      :: diameter
    dsl4jsb_Real2D_onChunk      :: height
    dsl4jsb_Real2D_onChunk      :: cohort_age
    dsl4jsb_Real2D_onChunk      :: mortality_rate
    dsl4jsb_Real2D_onChunk      :: delta_dens_ind
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: an_boc_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: net_growth_tvegdyn_mavg
    dsl4jsb_Real2D_onChunk      :: lai_tvegdyn_mavg
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
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_RAD_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    dsl4jsb_Get_mt2L2D(VEG_BGCM_LITTERFALL_ID, veg_litterfall_mt)
    dsl4jsb_Get_mt1L2D(VEG_BGCM_ESTABLISHMENT_ID, veg_establishment_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! Q_RAD_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, rfr_ratio_boc_tvegdyn_mavg) ! in
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)         ! in
    ! SPQ_ 2D
    dsl4jsb_Get_var2D_onChunk(SPQ_, w_soil_root_theta)          ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, do_cohort_harvest)          ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, dens_ind)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, diameter)                   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, height)                     ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, cohort_age)                 ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, mortality_rate)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, delta_dens_ind)             ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)            ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, an_boc_tvegdyn_mavg)        ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, net_growth_tvegdyn_mavg)    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, lai_tvegdyn_mavg)           ! in
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 Calculate foliage projective cover if vegetation dynamics are simulated, i.e., not needed for constant mortality
    !>
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
      CASE ("population","cohort")
        fpc(:) = calc_foliage_projective_cover_tile(nc, &
          &                                         lctlib%growthform, &
          &                                         lctlib%wood_density, &
          &                                         lctlib%k_latosa, &
          &                                         veg_pool_mt(ix_sap_wood, ixC, :), &
          &                                         dens_ind(:), &
          &                                         height(:), &
          &                                         diameter(:))
    ENDSELECT

    !>  1.1 check stand harvest year (if this functionality is enabled: model%config\%flag_stand_harvest)
    !>
    ! but enable harvesting only after transient spinup, i.e., when l_transient_spinup=FALSE
    IF (model%config%flag_stand_harvest .AND. .NOT. model%config%l_transient_spinup) THEN
      model%config%l_do_stand_replacing_harvest = .FALSE.                                      ! per default this flag is FALSE
      ! set stand_harvest to TRUE for the 1st timestep of the stand_replacing_year only
      IF (forcing_options_quincy%forcing_year .EQ. model%config%stand_replacing_year .AND. &
        & forcing_options_quincy%forcing_doy  .EQ. 1                                 .AND. &
        & forcing_options_quincy%forcing_hour .LT. eps4                                     ) THEN
        model%config%l_do_stand_replacing_harvest = .TRUE.
      END IF
    END IF

    !>  1.2 Increase age of the cohort [years] & check if cohort harvest is reached
    !>
    ! apply this only after spinup of transient mode has finished
    IF (.NOT. model%config%l_transient_spinup) THEN
      SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
        CASE ("cohort")
          cohort_age(:)        = cohort_age(:) + model%config%dtime / (one_year * one_day)
          do_cohort_harvest(:) = false                                               ! per default it is set to FALSE
          WHERE (cohort_age(:) > dsl4jsb_Config(VEG_)%cohort_harvest_interval)
            do_cohort_harvest(:) = true
            cohort_age(:)        = 0.0_wp
          ENDWHERE
      ENDSELECT
    ENDIF

    !>2.0 Calculate mortality rate and veg\%litterfall
    !>
    ! retruns constant mortality_rate if veg_dynamics_scheme == none \n
    ! but calcs mortality_rate dynamically if veg_dynamics_scheme == population or cohort
    CALL calc_veg_mortality( &
      & nc                                              , & ! in
      & model%config%elements_index_map(:)              , &
      & model%config%is_element_used(:)                 , &
      & lctlib%growthform                               , &
      & lctlib%k1_mort_greff                            , &
      & TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme)  , &
      & an_boc_tvegdyn_mavg(:)                          , &
      & net_growth_tvegdyn_mavg(:)                      , &
      & lai_tvegdyn_mavg(:)                             , &
      & fpc(:)                                          , &
      & veg_pool_mt(:,:,:)                              , & ! in
      & veg_litterfall_mt(:,:,:)                        , & ! inout
      & mortality_rate(:)                               )   ! out

    !>  2.1 Calculate stand harvest_rate and veg\%litterfall
    !>
    ! if either logical is true
    ! @TODO update_veg_dynamics() do_cohort_harvest dimension should be used as nc
    IF (model%config%l_do_stand_replacing_harvest .OR. do_cohort_harvest(ics) > test_false_true) THEN
      CALL calc_veg_harvest( &
        & nc, &                                 ! in
        & model%config%elements_index_map(:), &
        & model%config%is_element_used(:), &
        & veg_pool_mt(:,:,:), &                 ! in
        & veg_litterfall_mt(:,:,:), &           ! inout
        & mortality_rate(:) )                   ! inout
    ENDIF

    !>3.0 Calculate establishment_rate and veg\%establishment
    !>
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
    !CASE ("cohort")
      ! no establishment in cohort mode, but after cohort harvesting \n
      !  the dens_ind is set to 1.0 and veg_pool_mt(ix_reserve gets some amount of CNP same as during the veg_init \n
      !  in the routine update_veg_pools section 3.0
    CASE ("population")
      IF(.NOT. model%config%l_do_stand_replacing_harvest) THEN
        CALL calc_veg_establishment( &
          & nc                                  , & ! in
          & model%config%elements_index_map(:)  , &
          & model%config%is_element_used(:)     , &
          & lctlib%tau_seed_est                 , &
          & w_soil_root_theta(:)                , &
          & growing_season(:)                   , &
          & t_air_week_mavg(:)                  , &
          & fpc(:)                              , &
          & rfr_ratio_boc_tvegdyn_mavg(:)       , &
          & veg_pool_mt(ix_seed_bed,:,:)        , & ! in
          & veg_establishment_mt(:,:)           )   ! inout
      ENDIF
    ENDSELECT

    !>4.0 Calculate change in individual density from mortality and establishment
    !>
    ! only if no stand/cohort harvesting is done
    SELECT CASE (TRIM(dsl4jsb_Config(VEG_)%veg_dynamics_scheme))
    CASE ("population","cohort")
      !
      ! @TODO update_veg_dynamics() do_cohort_harvest dimension should be used as nc
      !
      IF(.NOT. model%config%l_do_stand_replacing_harvest .AND. do_cohort_harvest(ics) < test_false_true) THEN
        IF(lctlib%growthform == ITREE) THEN
          delta_dens_ind(:) = calc_delta_dens_ind(lctlib%seed_size            , &
            &                                     dens_ind(:)                 , &
            &                                     mortality_rate(:)           , &
            &                                     veg_establishment_mt(ixC,:) )
        ELSE
          delta_dens_ind(:) = 0.0_wp
        ENDIF
      ENDIF
    ENDSELECT

  END SUBROUTINE update_veg_dynamics

  ! ======================================================================================================= !
  !>calculates current foliage projective cover for a tile
  !>
  !>  follows Sitch et al. 2003, eq 8
  !>
  FUNCTION calc_foliage_projective_cover_tile( &
    & nc                      , &
    & lctlib_growthform       , &
    & lctlib_wood_density     , &
    & lctlib_k_latosa         , &
    & veg_pool_sap_wood_carbon, &
    & dens_ind                , &
    & height                  , &
    & diameter)                 &
    & RESULT(fpc)

    USE mo_veg_constants,         ONLY: max_crown_area, min_diameter, k_crown_area, k_rp, k_fpc, ITREE
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in) :: nc                          !< dimensions
    INTEGER,                  INTENT(in) :: lctlib_growthform           !< lctlib paramter
    REAL(wp),                 INTENT(in) :: lctlib_wood_density         !< lctlib paramter
    REAL(wp),                 INTENT(in) :: lctlib_k_latosa             !< lctlib paramter
    REAL(wp), DIMENSION(nc),  INTENT(in) :: veg_pool_sap_wood_carbon, & !< sapwood mass (mol C / m2)
                                            dens_ind                , & !< tree density (#/m2)
                                            height                  , & !< tree height (m)
                                            diameter                    !< tree diameter (m)
    REAL(wp), DIMENSION(nc)              :: fpc                         !< foliage projective cover (fraction)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp), DIMENSION(nc)                     :: crown_area
    REAL(wp), DIMENSION(nc)                     :: lai_ind
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_foliage_projective_cover_tile'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 ...
    !>
    WHERE((lctlib_growthform == ITREE) .AND. (dens_ind(:) > eps4) .AND. (height(:) > eps4))
      ! crown area derived from diameter following Sitch et al. 2003; eq 4
      crown_area(:) = MIN(max_crown_area, k_crown_area * MAX(min_diameter,diameter(:)) ** k_rp)
      ! maximum annual, individual LAI, derived from the sapwood area of an individual tree, devided by
      ! its crown area
      lai_ind(:)    = veg_pool_sap_wood_carbon(:) / dens_ind(:) / lctlib_wood_density / height(:) * lctlib_k_latosa / &
        &             crown_area(:)
      ! foliage projective cover is then the product of density and individual crown area
      fpc(:)        = crown_area(:) * dens_ind(:) * (1._wp - exp(-k_fpc * lai_ind(:)))
    ELSEWHERE
      lai_ind(:) = 0.0_wp
      fpc(:)     = 0.0_wp
    ENDWHERE
  END FUNCTION calc_foliage_projective_cover_tile

  ! ======================================================================================================= !
  !>calculate vegetation mortality
  !>
  SUBROUTINE calc_veg_mortality( &
    & nc                        , &
    & elements_index_map        , &
    & is_element_used           , &
    & lctlib_growthform         , &
    & lctlib_k1_mort_greff      , &
    & veg_dynamics_scheme       , &
    & an_boc_tvegdyn_mavg       , &
    & net_growth_tvegdyn_mavg   , &
    & lai_tvegdyn_mavg          , &
    & fpc                       , &
    & veg_pool_mt               , &
    & veg_litterfall_mt         , &
    & mortality_rate )

    USE mo_veg_constants,         ONLY: min_greff, k2_mort_greff, k3_mort_greff, fpc_max, &
      &                                 ITREE, IGRASS, background_mort_rate_tree, background_mort_rate_grass
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in)    :: nc                         !< dimensions
    INTEGER,                  INTENT(in)    :: elements_index_map(:)      !< map bgcm element ID -> IDX
    LOGICAL,                  INTENT(in)    :: is_element_used(:)         !< is element in 'elements_index_map' used
    INTEGER,                  INTENT(in)    :: lctlib_growthform          !< lctlib paramter
    REAL(wp),                 INTENT(in)    :: lctlib_k1_mort_greff       !< lctlib paramter
    CHARACTER(len=*),         INTENT(in)    :: veg_dynamics_scheme        !< vegetation dynamics: none population cohort
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: an_boc_tvegdyn_mavg        !< long-term net C balance at the bottom of the canopy (micro-mol / m2 / s)
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: net_growth_tvegdyn_mavg    !< long-term net growth of the plant (micro-mol / m2 / s)
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: lai_tvegdyn_mavg           !< long-term LAI of the plant
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: fpc                        !< foliage projective cover
    REAL(wp),                 INTENT(in)    :: veg_pool_mt(:,:,:)         !< bgcm veg_pool
    REAL(wp),                 INTENT(inout) :: veg_litterfall_mt(:,:,:)   !< bgcm flux: veg_litterfall
    REAL(wp), DIMENSION(nc),  INTENT(out)   :: mortality_rate             !< current mortality rate (1/timestep)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ic                             !< loop over point of chunk
    INTEGER                     :: ielem                          !< loop over bgcm elements
    INTEGER                     :: ix_elem                        !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc)     :: boc_cbalance                   !< C balance at bottom of canopy
    REAL(wp), DIMENSION(nc)     :: greff                          !< long-term growth efficiency
    REAL(wp), DIMENSION(nc)     :: fscal                          !< ...
    REAL(wp), DIMENSION(nc)     :: greff_mortality_rate           !< ...
    REAL(wp), DIMENSION(nc)     :: self_thinning_mortality_rate   !< ...
    REAL(wp), DIMENSION(nc)     :: bg_mortality_rate              !< ... of grasses or trees
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_mortality'
    ! ----------------------------------------------------------------------------------------------------- !

    !> 0.9 init local var
    !>
    boc_cbalance(:)                 = 0.0_wp
    greff(:)                        = 0.0_wp
    fscal(:)                        = 0.0_wp
    greff_mortality_rate(:)         = 0.0_wp
    self_thinning_mortality_rate(:) = 0.0_wp
    bg_mortality_rate(:)            = 0.0_wp

    ! differ between veg_dynamics_schemes: population/cohort & none
    SELECT CASE (veg_dynamics_scheme)
    CASE ("population","cohort")
      !>1.0 mortality related to growth efficiency mortality
      !>
      DO ic = 1,nc
        IF (lai_tvegdyn_mavg(ic) > eps8) THEN
          ! growth efficiency per unit leaf area in mol C / m2 LAI / yr
          !   the application of 'MAX()' avoids negative values of greff_mortality_rate(ic) and mortality_rate(ic)
          greff(ic) = MAX(net_growth_tvegdyn_mavg(ic) / lai_tvegdyn_mavg(ic) * one_day * one_year / 1.e6_wp, 0.0_wp)
          !>  1.1 deduce mortality from asymptoptic mortality rate given growth efficiency,
          !>     below a minimum threshold mortality rises to 100%
          !>
          !     @TODO the k1_mort_greff of PFT 3 (TrBR) is modified (increased compared to other PFT)
          !             to reflect disturbance (C loss) due to dry periods
          !             otherwise TrBR trees would die because of too strong maintanaince respiration during dry season,
          !           this "static" value of k1_mort_greff per PFT may be replaced by future implementations
          !             of explicit disturbance regimes
          IF (greff(ic) > min_greff) THEN
            greff_mortality_rate(ic) = lctlib_k1_mort_greff / (k2_mort_greff * greff(ic) + 1.0_wp)
          ELSE
            greff_mortality_rate(ic) = 1.0_wp / (k3_mort_greff * greff(ic) + 1.0_wp)
          END IF
        END IF
      END DO
      !>  1.2 self-thinning related mortality
      !>
      self_thinning_mortality_rate(:) = MAX(fpc(:) - fpc_max, 0.0_wp)
      ! ! depends on the carbon balance of lowest canopy layer, i.e. net photosynthesis minus average construction cost \n
      ! ! in mol / m2 LAI / year \n
      ! ! transformed to a mortality estimate using a Weibull function
      ! boc_cbalance(:) = an_boc_tvegdyn_mavg(:) * one_day * one_year / 1.e6_wp &
      !   &               - ( 1._wp + fresp_growth ) * 1.0_wp/lctlib_tau_leaf/lctlib_sla
      ! WHERE(boc_cbalance(:) < (-eps4))
      !   mortality_rate(:) = mortality_rate(:) &
      !     &                 + ( 1._wp - exp ( - ( - lambda_mort_light * boc_cbalance(:) ) ** k_mort_light))
      ! END WHERE

      !>  1.3 quasi-stochastic background mortality
      !>
      SELECT CASE (lctlib_growthform)
      CASE (ITREE)
        bg_mortality_rate(:) = background_mort_rate_tree
      CASE (IGRASS)
        bg_mortality_rate(:) = background_mort_rate_grass
      ENDSELECT
      !>  1.4 total mortality is the sum of all mortality terms, and limited to one
      !>
      mortality_rate(:) = MIN(greff_mortality_rate(:) + self_thinning_mortality_rate(:) + bg_mortality_rate(:), 1.0_wp) &
        &                 / one_year / one_day * dtime
    CASE ("none")
      SELECT CASE (lctlib_growthform)
      CASE (ITREE)
        mortality_rate(:) = ((background_mort_rate_tree /  one_year) / one_day) * dtime
      CASE (IGRASS)
        mortality_rate(:) = ((background_mort_rate_grass / one_year) / one_day) * dtime
      ENDSELECT
    ENDSELECT ! veg_dynamics_scheme

    !>2.0 calculate litter fall from mortality
    !>
    ! @NOTE for all tissues BUT NOT seed_bed
    DO ic = 1,nc
      IF (veg_pool_mt(ix_leaf, ixC, ic) > eps8) THEN
        fscal(ic) = (veg_pool_mt(ix_leaf, ixC, ic) - veg_litterfall_mt(ix_leaf, ixC, ic)) / veg_pool_mt(ix_leaf, ixC, ic)
      ELSE
        fscal(ic) = 0.0_wp
      END IF
      IF (mortality_rate(ic) < fscal(ic)) THEN
        fscal(ic) = mortality_rate(ic)
      END IF
    END DO
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_litterfall_mt(ix_leaf, ix_elem, :)        = veg_litterfall_mt(ix_leaf, ix_elem, :) &
          &                                             + veg_pool_mt(ix_leaf, ix_elem, :) * fscal(:)
        veg_litterfall_mt(ix_fine_root, ix_elem, :)   = veg_litterfall_mt(ix_fine_root, ix_elem, :) &
          &                                             + veg_pool_mt(ix_fine_root, ix_elem, :) * mortality_rate(:)
        veg_litterfall_mt(ix_coarse_root, ix_elem, :) = veg_litterfall_mt(ix_coarse_root, ix_elem, :) &
          &                                             + veg_pool_mt(ix_coarse_root, ix_elem, :) * mortality_rate(:)
        veg_litterfall_mt(ix_sap_wood, ix_elem, :)    = veg_litterfall_mt(ix_sap_wood, ix_elem, :) &
          &                                             + veg_pool_mt(ix_sap_wood, ix_elem, :) * mortality_rate(:)
        veg_litterfall_mt(ix_heart_wood, ix_elem, :)  = veg_litterfall_mt(ix_heart_wood, ix_elem, :) &
          &                                             + veg_pool_mt(ix_heart_wood, ix_elem, :) * mortality_rate(:)
        veg_litterfall_mt(ix_labile, ix_elem, :)      = veg_litterfall_mt(ix_labile, ix_elem, :) &
          &                                             + veg_pool_mt(ix_labile, ix_elem, :) * mortality_rate(:)
        veg_litterfall_mt(ix_reserve, ix_elem, :)     = veg_litterfall_mt(ix_reserve, ix_elem, :) &
          &                                             + veg_pool_mt(ix_reserve, ix_elem, :) * mortality_rate(:)
        veg_litterfall_mt(ix_fruit, ix_elem, :)       = veg_litterfall_mt(ix_fruit, ix_elem, :) &
          &                                             + veg_pool_mt(ix_fruit, ix_elem, :) * mortality_rate(:)
      END IF
    END DO
  END SUBROUTINE calc_veg_mortality

  ! ======================================================================================================= !
  !>calculate vegetation harvest
  !>
  SUBROUTINE calc_veg_harvest( &
    & nc                        , &
    & elements_index_map        , &
    & is_element_used           , &
    & veg_pool_mt               , &
    & veg_litterfall_mt         , &
    & mortality_rate)

    USE mo_veg_constants,         ONLY: wood_extraction_rate
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in)    :: nc                       !< dimensions
    INTEGER,                  INTENT(in)    :: elements_index_map(:)    !< map bgcm element ID -> IDX
    LOGICAL,                  INTENT(in)    :: is_element_used(:)       !< is element in 'elements_index_map' used
    REAL(wp),                 INTENT(in)    :: veg_pool_mt(:,:,:)       !< vegetation pool
    REAL(wp),                 INTENT(inout) :: veg_litterfall_mt(:,:,:) !< vegetation litterfall
    REAL(wp), DIMENSION(nc),  INTENT(inout) :: mortality_rate           !< background_mort_rate_tree/grass converted to 1/timestep
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ielem             !< loop over bgcm elements
    INTEGER                     :: ix_elem           !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc)     :: rest_mortality_rate
    REAL(wp)                    :: export_rate
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_harvest'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 calculate rest mortality to achieve stand-replacing harvest
    !>
    rest_mortality_rate(:) = 1._wp - mortality_rate(:)
    ! ...
    ! loop over bgcm elements and soil layers
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_litterfall_mt(ix_leaf, ix_elem, :)        = veg_litterfall_mt(ix_leaf, ix_elem, :) &
          &                                             + veg_pool_mt(ix_leaf, ix_elem, :) &
          &                                             * rest_mortality_rate(:)
        veg_litterfall_mt(ix_fine_root, ix_elem, :)   = veg_litterfall_mt(ix_fine_root, ix_elem, :) &
          &                                             + veg_pool_mt(ix_fine_root, ix_elem, :) &
          &                                             * rest_mortality_rate(:)
        veg_litterfall_mt(ix_coarse_root, ix_elem, :) = veg_litterfall_mt(ix_coarse_root, ix_elem, :) &
          &                                             + veg_pool_mt(ix_coarse_root, ix_elem, :) &
          &                                             * rest_mortality_rate(:)
        veg_litterfall_mt(ix_sap_wood, ix_elem, :)    = veg_litterfall_mt(ix_sap_wood, ix_elem, :) &
          &                                             + veg_pool_mt(ix_sap_wood, ix_elem, :) &
          &                                             * (1._wp - wood_extraction_rate) &
          &                                             * rest_mortality_rate(:)
        veg_litterfall_mt(ix_heart_wood, ix_elem, :)  = veg_litterfall_mt(ix_heart_wood, ix_elem, :) &
          &                                             + veg_pool_mt(ix_heart_wood, ix_elem, :) &
          &                                             * (1._wp - wood_extraction_rate) &
          &                                             * rest_mortality_rate(:)
        veg_litterfall_mt(ix_labile, ix_elem, :)      = veg_litterfall_mt(ix_labile, ix_elem, :) &
          &                                             + veg_pool_mt(ix_labile, ix_elem, :) &
          &                                             * (1._wp - wood_extraction_rate) &
          &                                             * rest_mortality_rate(:)
        veg_litterfall_mt(ix_reserve, ix_elem, :)     = veg_litterfall_mt(ix_reserve, ix_elem, :) &
          &                                             + veg_pool_mt(ix_reserve, ix_elem, :) &
          &                                             * (1._wp - wood_extraction_rate) &
          &                                             * rest_mortality_rate(:)
        veg_litterfall_mt(ix_fruit, ix_elem, :)       = veg_litterfall_mt(ix_fruit, ix_elem, :) &
          &                                             + veg_pool_mt(ix_fruit, ix_elem, :) &
          &                                             * rest_mortality_rate(:)
      END IF
    END DO
    ! ...
    mortality_rate(:) = 1.0_wp
  END SUBROUTINE calc_veg_harvest

  ! ======================================================================================================= !
  !>calculate vegetation establishment
  !>
  SUBROUTINE calc_veg_establishment( &
    & nc                          , &
    & elements_index_map          , &
    & is_element_used             , &
    & lctlib_tau_seed_est         , &
    & w_soil_root_theta           , &
    & growing_season              , &
    & t_air_week_mavg             , &
    & fpc                         , &
    & rfr_ratio_boc_tvegdyn_mavg  , &
    & veg_pool_mt_seed_bed        , &
    & veg_establishment_mt         )

    USE mo_q_rad_parameters,            ONLY: rfr_ratio_toc
    USE mo_veg_constants,               ONLY: fpc_max, lambda_est_temp, k_est_temp, lambda_est_moist, k_est_moist
    USE mo_jsb_physical_constants,      ONLY: Tzero
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,                  INTENT(in)    :: nc                           !< dimensions
    INTEGER,                  INTENT(in)    :: elements_index_map(:)        !< map bgcm element ID -> IDX
    LOGICAL,                  INTENT(in)    :: is_element_used(:)           !< is element in 'elements_index_map' used
    REAL(wp),                 INTENT(in)    :: lctlib_tau_seed_est          !< lctlib paramter
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: w_soil_root_theta            !< plant available water [fraction of maximum]
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: growing_season               !< growing season
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: t_air_week_mavg              !< weekly air temperature [K]
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: fpc                          !< folage projective cover
    REAL(wp), DIMENSION(nc),  INTENT(in)    :: rfr_ratio_boc_tvegdyn_mavg   !< red-farred ratio at the bottom of the canopy
    REAL(wp),                 INTENT(in)    :: veg_pool_mt_seed_bed(:,:)    !< bgcm veg_pool: seed_bed in vegetation pool
    REAL(wp),                 INTENT(inout) :: veg_establishment_mt(:,:)    !< bgcm flux: vegetation establishment
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                     :: ielem              !< loop over bgcm elements
    INTEGER                     :: ix_elem            !< index of element in bgcm, used for looping
    REAL(wp), DIMENSION(nc)     :: establishment_rate
    REAL(wp), DIMENSION(nc)     :: tc
    REAL(wp), DIMENSION(nc)     :: rfr_ratio
    REAL(wp), DIMENSION(nc)     :: flim_light
    REAL(wp), DIMENSION(nc)     :: flim_moist
    REAL(wp), DIMENSION(nc)     :: flim_temp
    INTEGER                     :: icanopy
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_veg_establishment'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 establishment rate
    !>

    !>  1.1 light limitation
    !>
    flim_light(:) = MAX(fpc_max - fpc(:), 0.0_wp)

    ! depending on the red to far-red ratio of the light at the bottom of the canopy
    ! flim_light(:) = EXP(-(lctlib_lambda_est_light * &
    !   &             (rfr_ratio_toc-rfr_ratio_boc_tvegdyn_mavg(:)))**lctlib_k_est_light)

    !> 1.2 temperature and soil moisture limitation
    !!
    !! @NOTE  OBS: should be top soil moisture!
    WHERE((t_air_week_mavg(:) - Tzero) < 0.0_wp .OR. growing_season(:) < test_false_true)
      flim_temp(:) = 0.0_wp
    ELSEWHERE
      flim_temp(:) = 1.0_wp - EXP(-(lambda_est_temp * (t_air_week_mavg(:) - Tzero)) ** k_est_temp)
    END WHERE
    flim_moist(:)  = 1.0_wp - EXP(-(lambda_est_moist * w_soil_root_theta(:)) ** k_est_moist)

    !>  1.3 actual establishment rate and associated matter flux from seed bed to reserve pool
    !>
    establishment_rate(:) = flim_light(:) * flim_temp(:) * flim_moist(:) &
      &                     * 1.0_wp / lctlib_tau_seed_est / one_day / one_year * dtime
    WHERE(establishment_rate(:) < eps8)
      establishment_rate(:) = 0.0_wp
    ENDWHERE
    ! calc vegetation establishment
    ! loop over bgcm elements
    DO ielem = FIRST_ELEM_ID, LAST_ELEM_ID
      IF (is_element_used(ielem)) THEN
        ix_elem = elements_index_map(ielem)    ! get element index in bgcm
        veg_establishment_mt(ix_elem, :) = veg_pool_mt_seed_bed(ix_elem, :) * establishment_rate(:)
      END IF
    END DO
  END SUBROUTINE calc_veg_establishment

  ! ======================================================================================================= !
  !>calculate change in tree density
  !>
  PURE ELEMENTAL FUNCTION calc_delta_dens_ind( &
    & lctlib_seed_size         , &
    & dens_ind                 , &
    & mortality_rate           , &
    & veg_establishment_carbon)  &
    & RESULT (delta_dens_ind)

    REAL(wp), INTENT(in) :: lctlib_seed_size            !< seed size to convert C into individuals, lctlib parameter
    REAL(wp), INTENT(in) :: dens_ind                    !< current individuum density (#/m2)
    REAL(wp), INTENT(in) :: mortality_rate              !< mortality per timestep
    REAL(wp), INTENT(in) :: veg_establishment_carbon    !< vegetation establishment measured in C
    REAL(wp)             :: delta_dens_ind              !< change in individual density
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_delta_density_individuals'
    ! ----------------------------------------------------------------------------------------------------- !

    delta_dens_ind = (veg_establishment_carbon / lctlib_seed_size) - (mortality_rate * dens_ind)
  END FUNCTION calc_delta_dens_ind

#endif
END MODULE mo_q_veg_dynamics
