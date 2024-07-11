!> QUINCY calculate vegetation phenology
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
!>#### determine onset and end of the growing season across plant functional types (PFT)
!>
MODULE mo_q_pheno_process
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: true, false, test_false_true
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_q_phenology

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_pheno_process'

CONTAINS

  ! ======================================================================================================= !
  !> update phenology - QUINCY
  !>   determine onset and end of the growing season for some PFT
  !>
  !>   PFT: evergreen, raingreen, cold-deciduous trees and perennial grasses
  !>
  !>  @TODO
  !>  OBS: currently uncalibrated, in particular rain green phenology
  !>  Dev: Have more cost-benefit in decision to start/stop growing season
  !>      what about seasonality beyond off/on behaviour? -> implies some collaboration between
  !>      turnover rates and phenology -> should the two be linked?
  SUBROUTINE update_q_phenology(tile, options)
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: A2L_, Q_PHENO_, Q_ASSIMI_, VEG_
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_PHENO_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER   :: model                  !< the model
    TYPE(t_lctlib_element), POINTER   :: lctlib                 !< land-cover-type library - parameter across pft's
    INTEGER                           :: iblk, ics, ice, nc     !< dimensions
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_q_phenology'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Real2D_onChunk      :: t_air
    ! Q_PHENO_ 2D
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: gdd
    dsl4jsb_Real2D_onChunk      :: nd_dormance
    ! Q_ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    ! VEG_ 2D
    dsl4jsb_Real2D_onChunk      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: t_soil_srf_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_PHENO_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN !< run task only if the present tile is a pft
    ! ----------------------------------------------------------------------------------------------------- !
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (lctlib%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_ 2D
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                            ! in
    ! Q_PHENO_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)               ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, gdd)                          ! inout
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, nd_dormance)                  ! inout
    ! Q_ASSIMI_ 2D
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tphen_mavg)     ! in
    ! VEG_ 2D
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_soil_srf_tphen_mavg)            ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                    ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)        ! inout
    ! ----------------------------------------------------------------------------------------------------- !

    !> 1.0 Check for start and end of growing season, and update both gdd and nd_dormance counter
    !>
    CALL calc_phenology(                     &
      & lctlib%phenology_type              , &  ! lctlib
      & lctlib%beta_soil_flush             , &
      & lctlib%cn_leaf                     , &
      & lctlib%np_leaf                     , &
      & lctlib%gdd_req_max                 , &
      & lctlib%k_gdd_dormance              , &
      & lctlib%beta_soil_senescence        , &
      & lctlib%t_air_senescence            , &
      & lctlib%min_leaf_age                , &  ! lctlib
      & t_air(:)                           , &  ! in
      & beta_soil_gs_tphen_mavg(:)         , &
      & t_air_tphen_mavg(:)                , &
      & t_soil_srf_tphen_mavg(:)           , &
      & t_air_week_mavg(:)                 , &
      & t_air_month_mavg(:)                , &
      & mean_leaf_age(:)                   , &
      & gpp_tlabile_mavg(:)                , &
      & maint_respiration_tlabile_mavg(:)  , &  ! in
      & growing_season(:)                  , &  ! inout
      & growth_req_n_tlabile_mavg(:)       , &
      & growth_req_p_tlabile_mavg(:)       , &
      & nd_dormance(:)                     , &
      & gdd(:)                             )    ! inout
  END SUBROUTINE update_q_phenology

  ! ======================================================================================================= !
  !> calculate phenology
  !>
  ELEMENTAL SUBROUTINE calc_phenology(  &
    & phenology_type                  , &
    & beta_soil_flush                 , &
    & cn_leaf                         , &
    & np_leaf                         , &
    & gdd_req_max                     , &
    & k_gdd_dormance                  , &
    & beta_soil_senescence            , &
    & t_air_senescence                , &
    & min_leaf_age                    , &
    & t_air                           , &
    & beta_soil_gs_tphen_mavg         , &
    & t_air_tphen_mavg                , &
    & t_soil_srf_tphen_mavg           , &
    & t_air_week_mavg                 , &
    & t_air_month_mavg                , &
    & mean_leaf_age                   , &
    & gpp_tlabile_mavg                , &
    & maint_respiration_tlabile_mavg  , &
    & growing_season                  , &
    & growth_req_n_tlabile_mavg       , &
    & growth_req_p_tlabile_mavg       , &
    & nd_dormance                     , &
    & gdd )
    ! ----------------------------------------------------------------------------------------------------- !
    USE mo_jsb_math_constants,      ONLY: dtime, one_day
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_jsb_impl_constants,      ONLY: true, false, test_false_true
    USE mo_q_pheno_constants,       ONLY: ievergreen, iraingreen, isummergreen, iperennial
    USE mo_q_pheno_parameters,      ONLY: gdd_t_air_threshold
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER,      INTENT(in)    :: phenology_type                   !< lctlib parameter
    REAL(wp),     INTENT(in)    :: beta_soil_flush                  !< lctlib parameter
    REAL(wp),     INTENT(in)    :: cn_leaf                          !< lctlib parameter
    REAL(wp),     INTENT(in)    :: np_leaf                          !< lctlib parameter
    REAL(wp),     INTENT(in)    :: gdd_req_max                      !< lctlib parameter
    REAL(wp),     INTENT(in)    :: k_gdd_dormance                   !< lctlib parameter
    REAL(wp),     INTENT(in)    :: beta_soil_senescence             !< lctlib parameter
    REAL(wp),     INTENT(in)    :: t_air_senescence                 !< lctlib parameter
    REAL(wp),     INTENT(in)    :: min_leaf_age                     !< lctlib parameter
    REAL(wp),     INTENT(in)    :: t_air                            !< air temp
    REAL(wp),     INTENT(in)    :: beta_soil_gs_tphen_mavg          !< soil ..
    REAL(wp),     INTENT(in)    :: t_air_tphen_mavg                 !< moving average air temperature over phenology-relevant period
    REAL(wp),     INTENT(in)    :: t_soil_srf_tphen_mavg            !< moving average soil surface temperature
    REAL(wp),     INTENT(in)    :: t_air_week_mavg                  !< moving average air temperature over a week
    REAL(wp),     INTENT(in)    :: t_air_month_mavg                 !< moving average air temperature over a month
    REAL(wp),     INTENT(in)    :: mean_leaf_age                    !< mean leaf age
    REAL(wp),     INTENT(in)    :: gpp_tlabile_mavg                 !< moving average gpp over tlabile
    REAL(wp),     INTENT(in)    :: maint_respiration_tlabile_mavg   !< moving average maintenance respiration over tlabile
    REAL(wp),     INTENT(inout) :: growing_season                   !< growing season Yes/No
    REAL(wp),     INTENT(inout) :: growth_req_n_tlabile_mavg        !< moles N required for a unit of C growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: growth_req_p_tlabile_mavg        !< moles P required for a unit of N growth to determine labile pool size
    REAL(wp),     INTENT(inout) :: nd_dormance                      !< number of days in dormancy since last growing season
    REAL(wp),     INTENT(inout) :: gdd                              !< growing degree days
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp)                    :: hlp1                             !< helper variable for tmp values
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_phenology'
    ! ----------------------------------------------------------------------------------------------------- !

    !>1.0 check growing season has started or ended
    !>
    ! If not in growing season, check whether growing season has started
    IF(growing_season < test_false_true) THEN
      SELECT CASE (phenology_type)
      CASE (ievergreen)
        ! nothing to do here
      CASE (iraingreen)
        ! if soil moisture is above specified threshold
        IF(beta_soil_gs_tphen_mavg > beta_soil_flush) THEN
          growing_season            = true
          growth_req_n_tlabile_mavg = 1._wp / cn_leaf
          growth_req_p_tlabile_mavg = 1._wp / np_leaf
        ENDIF
      CASE (isummergreen)
        ! if GDD requirement is exceeded. The latter depends on the length of the
        !  dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
        hlp1 = gdd_req_max * EXP( - k_gdd_dormance * nd_dormance )
        IF(gdd > hlp1) THEN
           growing_season            = true
           growth_req_n_tlabile_mavg = 1._wp / cn_leaf
           growth_req_p_tlabile_mavg = 1._wp / np_leaf
        ENDIF
      CASE (iperennial)
        ! if temperature and soil moisture are not limiting growth
        hlp1 = gdd_req_max * EXP( - k_gdd_dormance * nd_dormance )
        IF(gdd > hlp1 .AND. beta_soil_gs_tphen_mavg > beta_soil_flush) THEN
          growing_season            = true
          growth_req_n_tlabile_mavg = 1._wp / cn_leaf
          growth_req_p_tlabile_mavg = 1._wp / np_leaf
        ENDIF
      END SELECT
    ! If in the growing season, check whether growing season has ended
    ELSE
      SELECT CASE (phenology_type)
      CASE (ievergreen)
        ! nothing to do here
      CASE (iraingreen)
        ! if soil moisture is below given threshold AND leaves have been out a given number of days
        IF(beta_soil_gs_tphen_mavg < beta_soil_senescence .AND. mean_leaf_age > min_leaf_age ) THEN
          growing_season = false
        ENDIF
      CASE (isummergreen)
        ! if temperature are declining below a certain threshold AND
        !  leaves have been out a given number of days
        IF(t_air_week_mavg < t_air_month_mavg .AND. t_air_tphen_mavg < t_air_senescence .AND. &
          & mean_leaf_age > min_leaf_age) THEN
          growing_season = false
        ENDIF
      CASE (iperennial)
        ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND
        !  leaves have been out a given number of days
        IF((t_air_tphen_mavg < t_air_senescence .OR. beta_soil_gs_tphen_mavg < beta_soil_senescence .OR. &
          & gpp_tlabile_mavg < maint_respiration_tlabile_mavg).AND. mean_leaf_age > min_leaf_age) THEN
          growing_season = false
        ENDIF
      END SELECT
    ENDIF

    !>2.0 modify phenology counters
    !>

    !>  2.1 gdd counter
    !>
    IF(t_air > gdd_t_air_threshold) THEN
      gdd = gdd + ((t_air - gdd_t_air_threshold) * dtime / one_day)
    ENDIF
    ! vegetation growth doesn't start when there is snow coverage and soil is frozen
    IF(t_air_tphen_mavg < Tzero .AND. t_soil_srf_tphen_mavg <= Tzero) THEN
      gdd = 0.0_wp
    ENDIF

    !>  2.2 nd_dormance counter
    !>
    IF (growing_season > test_false_true) THEN
      gdd         = 0.0_wp
      nd_dormance = 0.0_wp
    ELSE
      nd_dormance = nd_dormance + dtime / one_day
    ENDIF
  END SUBROUTINE calc_phenology

#endif
END MODULE mo_q_pheno_process
