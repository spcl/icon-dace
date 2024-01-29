!> module to determine the phenological state of a plant 
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_pheno_update_phenology
#ifndef __NO_JSBACH__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_impl_constants,  ONLY: true, false, test_false_true
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: update_phenology_quincy

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pheno_update_phenology'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  ! Main Task
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> determine onset and end of the growing season for some PFT
  !! 
  !! PFT: evergreen, raingreen, cold-deciduous trees and perennial grasses \n
  !! 
  !! @todo OBS: currently uncalibrated, in particular rain green phenology \n
  !!  Dev: Have more cost-benefit in decision to start/stop growing season \n
  !!      what about seasonality beyond off/on behaviour? -> implies some collaboration between
  !!      turnover rates and phenology -> should the two be linked? 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_phenology_quincy(tile, options)

    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_lctlib_class,        ONLY: t_lctlib_element
    USE mo_jsb_process_class,       ONLY: A2L_, ASSIMI_, PHENO_, RAD_ , VEG_
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_math_constants,      ONLY: dtime, one_day, one_year, eps4, eps8
    USE mo_veg_constants,           ONLY: k_leafon_canopy, max_leaf_shedding_rate
    USE mo_pheno_constants,         ONLY: ievergreen
    USE mo_veg_canopy,              ONLY: calc_canopy_layers

    ! Use of process configurations (t_PROC_config)
    !

    ! Use of process memories
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(PHENO_)
    dsl4jsb_Use_memory(ASSIMI_)
    dsl4jsb_Use_memory(RAD_)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER       :: model
    TYPE(t_lctlib_element), POINTER       :: lctlib       !< land-cover-type library - parameter across pft's
    TYPE(t_jsb_vgrid),      POINTER       :: vgrid_canopy ! Vertical grid
    INTEGER                               :: ncanopy      ! number of canopy layers
    INTEGER                               :: iblk, ics, ice, nc
      !! temporary solution for optional arguments for calc_canopy_layers
      REAL(wp), DIMENSION(options%nc) :: &
                                          t_air_tfrac_mavg                            , &
                                          t_air_tacclim_mavg                          , &
                                          press_srf_tfrac_mavg                        , &
                                          co2_mixing_ratio_tfrac_mavg                 , &
                                          ga_tfrac_mavg                               , &
                                          beta_air_tfrac_mavg                         , &
                                          beta_soa_tphen_mavg                         , &
                                          beta_soil_ps_tfrac_mavg                     , &
                                          beta_sinklim_ps_tfrac_mavg                  , &
                                          beta_soil_gs_tfrac_mavg
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':update_phenology_quincy'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    ! Declare pointers to variables in memory
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    ! PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    dsl4jsb_Real2D_onChunk      :: gdd
    dsl4jsb_Real2D_onChunk      :: nd_dormance
    dsl4jsb_Real2D_onChunk      :: lai
    dsl4jsb_Real2D_onChunk      :: lai_max                    ! site level specific max lai value, read from lai_clim input var
    ! ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    ! RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    ! VEG_
    dsl4jsb_Real2D_onChunk      :: t_air_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_week_mavg
    dsl4jsb_Real2D_onChunk      :: t_air_month_mavg
    dsl4jsb_Real2D_onChunk      :: mean_leaf_age
    dsl4jsb_Real2D_onChunk      :: gpp_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: maint_respiration_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_n_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: growth_req_p_tlabile_mavg
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_carbon
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_nitrogen
    dsl4jsb_Real2D_onChunk      :: veg_pool_leaf_phosphorus
    dsl4jsb_Real2D_onChunk      :: veg_growth_leaf_carbon
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_leaf_carbon
    dsl4jsb_Real3D_onChunk      :: lai_cl
    dsl4jsb_Real3D_onChunk      :: cumm_lai_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: fleaf_sunlit_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: leaf_nitrogen_cl
    dsl4jsb_Real3D_onChunk      :: fn_chl_cl
    dsl4jsb_Real3D_onChunk      :: fn_et_cl
    dsl4jsb_Real3D_onChunk      :: fn_rub_cl
    dsl4jsb_Real3D_onChunk      :: fn_pepc_cl
    dsl4jsb_Real3D_onChunk      :: fn_oth_cl
    dsl4jsb_Real2D_onChunk      :: t_jmax_opt_mavg
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (.NOT. tile%Is_process_active(PHENO_)) RETURN
    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ---------------------------
    ! 0.5 Get Memory
    ! Get process variables (Set pointers to variables in memory)
    model         => Get_model(tile%owner_model_id)
    lctlib        => model%lctlib(tile%lcts(1)%lib_id)
    vgrid_canopy  => Get_vgrid('canopy_layer_q_')
    ncanopy       =  vgrid_canopy%n_levels
    ! Get process config
    !
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(PHENO_, growing_season)                 ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_, gdd)                            ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_, nd_dormance)                    ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai)                            ! inout
    dsl4jsb_Get_var2D_onChunk(PHENO_, lai_max)                        ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                            ! in
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_, beta_soil_gs_tphen_mavg)       ! in
    ! ---------------------------
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_cl)                      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_cl)                      ! in
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_sunlit_tfrac_mavg_cl)           ! inout
    dsl4jsb_Get_var3D_onChunk(RAD_, ppfd_shaded_tfrac_mavg_cl)           ! inout
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_tphen_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_week_mavg)                  ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, t_air_month_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, mean_leaf_age)                    ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, gpp_tlabile_mavg)                 ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_tlabile_mavg)   ! in
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_n_tlabile_mavg)        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, growth_req_p_tlabile_mavg)        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_carbon)             ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_nitrogen)           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_pool_leaf_phosphorus)         ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_leaf_carbon)           ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_leaf_carbon)       ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, lai_cl)                           ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, cumm_lai_cl)                      ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_cl)                  ! in
    dsl4jsb_Get_var3D_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl)       ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, leaf_nitrogen_cl)                 ! out
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_chl_cl)                        ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_et_cl)                         ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_rub_cl)                        ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_pepc_cl)                       ! inout
    dsl4jsb_Get_var3D_onChunk(VEG_, fn_oth_cl)                        ! inout
    dsl4jsb_Get_var2D_onChunk(VEG_, t_jmax_opt_mavg)                  !   optional input for calc_canopy_layers



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 0.8 pass reasonable values to mavg variables (temporary solution)
    !! 
    ppfd_sunlit_tfrac_mavg_cl(:,:)  = ppfd_sunlit_cl(:,:) 
    ppfd_shaded_tfrac_mavg_cl(:,:)  = ppfd_shaded_cl(:,:) 
    fleaf_sunlit_tfrac_mavg_cl(:,:) = fleaf_sunlit_cl(:,:)

    !> 0.9 optional arguments of calc_canopy_layers (local variables)
    !! 
    t_air_tfrac_mavg(:)            = 0.0_wp
    t_air_tacclim_mavg(:)          = 0.0_wp
    press_srf_tfrac_mavg(:)        = 0.0_wp
    co2_mixing_ratio_tfrac_mavg(:) = 0.0_wp
    ga_tfrac_mavg(:)               = 0.0_wp
    beta_air_tfrac_mavg(:)         = 0.0_wp
    beta_soa_tphen_mavg(:)         = 0.0_wp
    beta_soil_ps_tfrac_mavg(:)     = 0.0_wp
    beta_sinklim_ps_tfrac_mavg(:)  = 0.0_wp
    beta_soil_gs_tfrac_mavg(:)     = 0.0_wp


    !> 1.0 Check for start and end of growing season, and update both gdd and nd_dormance counter
    !! 
    CALL calc_phenology(lctlib%phenology_type              , &  ! lctlib
                        lctlib%beta_soil_flush             , &
                        lctlib%cn_leaf                     , &
                        lctlib%np_leaf                     , &
                        lctlib%gdd_req_max                 , &
                        lctlib%k_gdd_dormance              , &
                        lctlib%beta_soil_senescence        , &
                        lctlib%t_air_senescence            , &
                        lctlib%min_leaf_age                , &  ! lctlib
                        t_air(:)                           , &  ! in
                        beta_soil_gs_tphen_mavg(:)         , &
                        t_air_tphen_mavg(:)                , &
                        t_air_week_mavg(:)                 , &
                        t_air_month_mavg(:)                , &
                        mean_leaf_age(:)                   , &
                        gpp_tlabile_mavg(:)                , &
                        maint_respiration_tlabile_mavg(:)  , &
                        growing_season(:)                  , &  ! inout
                        growth_req_n_tlabile_mavg(:)       , &
                        growth_req_p_tlabile_mavg(:)       , &
                        nd_dormance(:)                     , &
                        gdd(:)                             )    ! inout


    !> 2.0 QUINCY canopy model: calculate lai(:) using growing season and leaf_carbon
    !! 
    !! this code is a simplified version of the 'veg_update_pools CASE ('canopy')' code in QUINCY
    
    !> 2.1 leaf carbon growth/litterfall flux depending on growing season T/F
    !!
    !! use in calc_canopy_layers() the lai(:) calculated here in sec. 2.3
    ! where growing_season is true
    WHERE(growing_season(:) > test_false_true)
       veg_growth_leaf_carbon(:) = lai_max(:) / lctlib%sla * &
                                   (1._wp - exp( -k_leafon_canopy * &
                                    (lai_max(:) - lai(:)) / lai_max(:))) * &
                                   dtime / one_day 
    ELSEWHERE
       veg_litterfall_leaf_carbon(:) = MIN(veg_pool_leaf_carbon(:), & 
                                           lai(:) / lctlib%sla * max_leaf_shedding_rate * dtime / one_day * & 
                                           lai_max(:) / &
                                           MAX(eps4, lai(:)))
    ENDWHERE

    !> 2.2 calc leaf turnover for evergreen trees (needed, for e.g. mean_leaf_age, as evergreens are always within growing_season)
    !!
    IF(lctlib%phenology_type == ievergreen) THEN
      veg_litterfall_leaf_carbon(:) = veg_pool_leaf_carbon(:) / lctlib%tau_leaf / one_day / one_year * dtime
    ENDIF
    
    !> 2.3 update leaf C:N:P pools 
    !!
    veg_pool_leaf_carbon(:)     = veg_pool_leaf_carbon(:)   + veg_growth_leaf_carbon(:) - veg_litterfall_leaf_carbon(:) 
    veg_pool_leaf_nitrogen(:)   = veg_pool_leaf_carbon(:)   / lctlib%cn_leaf 
    veg_pool_leaf_phosphorus(:) = veg_pool_leaf_nitrogen(:) / lctlib%np_leaf 

    !> 2.4 calc lai
    !!
    lai(:)      = veg_pool_leaf_carbon(:) * lctlib%sla

    !> 2.5 calc mean_leaf_age (copied from QUINCY update_veg_pools)
    !!
    WHERE(veg_pool_leaf_carbon(:) > eps8)
        mean_leaf_age(:) = (mean_leaf_age(:) + dtime / one_day) * (1._wp - veg_growth_leaf_carbon(:) / veg_pool_leaf_carbon(:)) &
                           - mean_leaf_age(:) * veg_litterfall_leaf_carbon(:) / veg_pool_leaf_carbon(:)
      ELSEWHERE
        mean_leaf_age(:) = 0.0_wp
    END WHERE


    !> 3.0 canopy layers
    !! 
    !! copied from quincy:mo_veg_update_pools 
    !!
    !! quincy canopy model:
    !!      previously this routine was called in update_canopy_fluxes() but led to inconsistentcies between lai & lai_cl in RAD
    CALL calc_canopy_layers( nc                                             , & ! in
                             ncanopy                                        , &
                             vgrid_canopy%dz(:)                             , &
                             vgrid_canopy%lbounds(:)                        , &
                             vgrid_canopy%ubounds(:)                        , &
                             lctlib%ps_pathway                              , &
                             lctlib%k0_fn_struc                             , &
                             lctlib%fn_oth_min                              , &
                             lctlib%sla                                     , &
                             lctlib%np_leaf                                 , &
                             lctlib%gmin                                    , &
                             lctlib%g0                                      , &
                             lctlib%g1                                      , &
                             lctlib%t_jmax_omega                            , &
                             .FALSE. , & !dsl4jsb_Config(ASSIMI_)%flag_optimal_Nfraction , &
                             veg_pool_leaf_nitrogen(:) , &  !veg_pool%leaf%nitrogen(:)                      , &
                             lai(:)                                         , & ! in
                             ppfd_sunlit_tfrac_mavg_cl(:,:)                 , & ! inout
                             ppfd_shaded_tfrac_mavg_cl(:,:)                 , & ! inout
                             fleaf_sunlit_tfrac_mavg_cl(:,:)                , & ! inout
                             t_air_tfrac_mavg(:)                            , &   ! optional, declared local var & set to 0.0_wp
                             t_air_tacclim_mavg(:)                          , &   ! optional, declared local var & set to 0.0_wp
                             press_srf_tfrac_mavg(:)                        , &   ! optional, declared local var & set to 0.0_wp
                             co2_mixing_ratio_tfrac_mavg(:)                 , &   ! optional, declared local var & set to 0.0_wp
                             ga_tfrac_mavg(:)                               , &   ! optional, declared local var & set to 0.0_wp
                             beta_air_tfrac_mavg(:)                         , &   ! optional, declared local var & set to 0.0_wp
                             beta_soa_tphen_mavg(:)                         , &   ! optional, declared local var & set to 0.0_wp
                             beta_soil_ps_tfrac_mavg(:)                     , &   ! optional, declared local var & set to 0.0_wp
                             beta_sinklim_ps_tfrac_mavg(:)                  , &   ! optional, declared local var & set to 0.0_wp
                             beta_soil_gs_tfrac_mavg(:)                     , &   ! optional, declared local var & set to 0.0_wp
                             t_jmax_opt_mavg(:)                             , &   ! optional
                             lai_cl(:,:)                                    , & ! out
                             cumm_lai_cl(:,:)                               , & ! out
                             leaf_nitrogen_cl(:,:)                          , & ! out
                             fn_rub_cl(:,:)                                 , & ! inout
                             fn_et_cl(:,:)                                  , & ! inout
                             fn_pepc_cl(:,:)                                , & ! inout
                             fn_chl_cl(:,:)                                 , & ! inout
                             fn_oth_cl(:,:)                                 )   ! inout

  END SUBROUTINE update_phenology_quincy 


  !-----------------------------------------------------------------------------------------------------
  ! Sub Task of update_phenology_quincy
  ! 
  !-----------------------------------------------------------------------------------------------------
  !> calculate phenology 
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE calc_phenology(phenology_type                  , &  ! lctlib
                                      beta_soil_flush                 , &
                                      cn_leaf                         , &
                                      np_leaf                         , &
                                      gdd_req_max                     , &
                                      k_gdd_dormance                  , &
                                      beta_soil_senescence            , &
                                      t_air_senescence                , &
                                      min_leaf_age                    , &  ! lctlib
                                      t_air                           , &
                                      beta_soil_gs_tphen_mavg         , &
                                      t_air_tphen_mavg                , &
                                      t_air_week_mavg                 , &
                                      t_air_month_mavg                , &
                                      mean_leaf_age                   , &
                                      gpp_tlabile_mavg                , &
                                      maint_respiration_tlabile_mavg  , &
                                      growing_season                  , & ! inout
                                      growth_req_n_tlabile_mavg       , &
                                      growth_req_p_tlabile_mavg       , &
                                      nd_dormance                     , &
                                      gdd )

    USE mo_jsb_math_constants,      ONLY: dtime, one_day
    USE mo_jsb_physical_constants,  ONLY: tmelt  
    USE mo_pheno_constants,         ONLY: ievergreen, iraingreen, isummergreen, iperennial, &
                                          gdd_t_air_threshold

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                      INTENT(in)    :: phenology_type                      !< land-cover-type library parameter
    REAL(wp),                     INTENT(in)    :: beta_soil_flush                 , & !< land-cover-type library parameter
                                                   cn_leaf                         , & !< land-cover-type library parameter
                                                   np_leaf                         , & !< land-cover-type library parameter
                                                   gdd_req_max                     , & !< land-cover-type library parameter
                                                   k_gdd_dormance                  , & !< land-cover-type library parameter
                                                   beta_soil_senescence            , & !< land-cover-type library parameter
                                                   t_air_senescence                , & !< land-cover-type library parameter
                                                   min_leaf_age                        !< land-cover-type library parameter
    REAL(wp),                     INTENT(in)    :: t_air                           , & !< air temp
                                                   beta_soil_gs_tphen_mavg         , & !< soil ..
                                                   t_air_tphen_mavg                , & !< moving average air temperature over ...
                                                   t_air_week_mavg                 , & !< moving average air temperature over a week
                                                   t_air_month_mavg                , & !< moving average air temperature over a month
                                                   mean_leaf_age                   , & !< mean leaf age
                                                   gpp_tlabile_mavg                , & !< moving average gpp over tlabile
                                                   maint_respiration_tlabile_mavg      !< moving average maintenance respiration over tlabile
    REAL(wp),                     INTENT(inout) :: growing_season                      !< growing season Yes/No
    REAL(wp),                     INTENT(inout) :: growth_req_n_tlabile_mavg       , & !< moles N required for a unit of C growth to determine labile pool size
                                                   growth_req_p_tlabile_mavg       , & !< moles P required for a unit of N growth to determine labile pool size
                                                   nd_dormance                     , & !< number of days in dormancy since last growing season
                                                   gdd                                 !< growing degree days 
    ! ---------------------------
    ! 0.2 Local
    REAL(wp) :: hlp1                                      ! helper variable for tmp values
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':calc_phenology'



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !> 1.0 check growing season has started or ended
    !! 
    IF(growing_season < test_false_true) THEN ! If not in growing season, check whether growing season has started
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is above specified threshold
            IF(beta_soil_gs_tphen_mavg > beta_soil_flush) THEN 
                  growing_season = true
                  growth_req_n_tlabile_mavg = 1._wp / cn_leaf
                  growth_req_p_tlabile_mavg = 1._wp / np_leaf
            ENDIF
         CASE (isummergreen)     ! if GDD requirement is exceeded. The latter depends on the length of the 
                                 ! dormancy period (aka "chilling requirement"). Parameterisation follows ORCHIDEE
            hlp1 = gdd_req_max * exp( - k_gdd_dormance * nd_dormance ) 
            IF(gdd > hlp1) THEN 
                  growing_season = true
                  growth_req_n_tlabile_mavg = 1._wp / cn_leaf
                  growth_req_p_tlabile_mavg = 1._wp / np_leaf
            ENDIF
         CASE (iperennial)       ! if temperature and soil moisture are not limiting growth
            hlp1 = gdd_req_max * exp( - k_gdd_dormance * nd_dormance ) 
            IF(gdd > hlp1 .AND. &
               beta_soil_gs_tphen_mavg > beta_soil_flush)THEN 
                  growing_season = true
                  growth_req_n_tlabile_mavg = 1._wp / cn_leaf
                  growth_req_p_tlabile_mavg = 1._wp / np_leaf
            ENDIF
       END SELECT
    ELSE                         ! If in the growing season, check whether growing season has ended
       SELECT CASE (phenology_type)
         CASE (ievergreen)       ! nothing to do here
         CASE (iraingreen)       ! if soil moisture is below given threshold AND leaves have been out a given number of days
            IF(beta_soil_gs_tphen_mavg < beta_soil_senescence .AND. & 
               mean_leaf_age > min_leaf_age ) & 
                  growing_season = false
         CASE (isummergreen)     ! if temperature are declining below a certain threshold AND 
            ! leaves have been out a given number of days
            IF(t_air_week_mavg < t_air_month_mavg .AND. &
               t_air_tphen_mavg < t_air_senescence .AND. & 
               mean_leaf_age > min_leaf_age) & 
                  growing_season = false
         CASE (iperennial)       ! if carbon balance of the plant becomes negative OR leaves have been damaged by frost/drought AND 
                                 ! leaves have been out a given number of days
            IF((t_air_tphen_mavg < t_air_senescence .OR. & 
                beta_soil_gs_tphen_mavg < beta_soil_senescence .OR. &
                gpp_tlabile_mavg < maint_respiration_tlabile_mavg).AND. &
                mean_leaf_age > min_leaf_age) &
                  growing_season = false
       END SELECT
    ENDIF   

    !> 2.1 gdd counter
    !! 
    IF(t_air > gdd_t_air_threshold) THEN
       gdd = gdd + ((t_air - gdd_t_air_threshold) * dtime/one_day)
    ENDIF
    IF(t_air_tphen_mavg < tmelt) gdd = 0.0_wp

    !> 2.2 nd_dormance counter
    !! 
    IF (growing_season > test_false_true) THEN 
      gdd         = 0.0_wp
      nd_dormance = 0.0_wp
    ELSE
      nd_dormance = nd_dormance + dtime/one_day 
    ENDIF

  END SUBROUTINE calc_phenology 

#endif
END MODULE mo_pheno_update_phenology
