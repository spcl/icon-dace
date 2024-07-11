!> QUINCY radiation process interface
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
!>#### definition and init of tasks for the radiation process incl. update and aggregate routines
!>
MODULE mo_q_rad_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_exception,           ONLY: message, finish
  USE mo_jsb_model_class,     ONLY: t_jsb_model
  USE mo_jsb_class,           ONLY: Get_model
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,   ONLY: t_jsb_process
  USE mo_jsb_task_class,      ONLY: t_jsb_process_task, t_jsb_task_options

  USE mo_lnd_bgcm_idx
  USE mo_lnd_bgcm_store,          ONLY: t_lnd_bgcm_store
  USE mo_lnd_bgcm_store_class,    ONLY: VEG_BGCM_POOL_ID


  dsl4jsb_Use_processes Q_RAD_, VEG_, A2L_, Q_PHENO_
  dsl4jsb_Use_memory(Q_RAD_)
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(Q_PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_q_rad_tasks
#ifdef __QUINCY_STANDALONE__
  PUBLIC :: update_time_average_q_radiation
#endif

  ! ======================================================================================================= !
  !> q_radiation
  !>
  ! Type definition
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_q_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_radiation
    PROCEDURE, NOPASS :: Aggregate => aggregate_q_radiation
  END TYPE tsk_q_radiation
  ! Constructor interface
  INTERFACE tsk_q_radiation
    PROCEDURE Create_task_q_radiation
  END INTERFACE tsk_q_radiation

  ! ======================================================================================================= !
  !> time_average_radiation task
  !>
  ! TYPE definition
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_q_radiation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_q_radiation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_q_radiation !< Aggregates computed task variables
  END TYPE tsk_time_average_q_radiation
  ! Constructor interface
  INTERFACE tsk_time_average_q_radiation
    PROCEDURE Create_task_update_time_average_q_radiation                !< Constructor function for task
  END INTERFACE tsk_time_average_q_radiation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_rad_interface'

CONTAINS

  ! ======================================================================================================= !
  !> Constructors for tasks
  !>

  !> q_radiation
  !>
  FUNCTION Create_task_q_radiation(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_radiation::return_ptr)
    CALL return_ptr%Construct(name='q_radiation', process_id=Q_RAD_, owner_model_id=model_id)
  END FUNCTION Create_task_q_radiation

  !> time_average_q_radiation
  !>
  FUNCTION Create_task_update_time_average_q_radiation(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_q_radiation::return_ptr)
    CALL return_ptr%Construct(name='tavrg_q_radiation', process_id=Q_RAD_, owner_model_id=model_id)
  END FUNCTION Create_task_update_time_average_q_radiation

  ! ======================================================================================================= !
  !> Register tasks for q_radiation process
  !>
  SUBROUTINE Register_q_rad_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER, INTENT(in) :: model_id

    CALL this%Register_task(tsk_q_radiation(model_id))
    CALL this%Register_task(tsk_time_average_q_radiation(model_id))
  END SUBROUTINE Register_q_rad_tasks

  ! ======================================================================================================= !
  !>Implementation of "integrate": q_radiation task
  !>
  SUBROUTINE update_q_radiation(tile, options)
    USE mo_q_rad_process,   ONLY: q_radiation
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_q_radiation'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk = options%iblk
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_RAD_)) RETURN
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !

    CALL q_radiation(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_q_radiation

  ! ======================================================================================================= !
  !> aggregate: q_radiation task
  !>
  SUBROUTINE aggregate_q_radiation(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_RAD_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_q_radiation'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_RAD_)
    dsl4jsb_Get_memory(VEG_)

    ! Q_RAD_ 2D
    dsl4jsb_Aggregate_onChunk(Q_RAD_, sw_srf_net                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, swvis_srf_net             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, swnir_srf_net             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, rad_srf_net               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_vis                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_nir                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_vis_soil              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_nir_soil              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_vis_snow              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_nir_snow              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_vis_can               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, alb_nir_can               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_vis_soil             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_nir_soil             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_vis_can              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_nir_can              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, appfd                     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, rfr_ratio_boc             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, albedo                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, net_radiation             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, albedo_noon               , weighted_by_fract)
    ! Q_RAD_ 3D
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_sunlit_cl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_shaded_cl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_sunlit_vis_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_shaded_vis_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, fleaf_sunlit_vis_cl       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_sunlit_nir_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, arad_shaded_nir_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, fleaf_sunlit_nir_cl       , weighted_by_fract)
    ! VEG_ 3D
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_cl           , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE aggregate_q_radiation

  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for Q_RAD_
  !>
  SUBROUTINE update_time_average_q_radiation(tile, options)

    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8, eps1
    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_grid_class,          ONLY: t_jsb_vgrid
    USE mo_jsb_grid,                ONLY: Get_vgrid
    USE mo_lnd_time_averages        ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    dsl4jsb_Use_processes A2L_, VEG_, Q_RAD_, Q_PHENO_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(VEG_)
    dsl4jsb_Use_memory(Q_RAD_)
    dsl4jsb_Use_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_vgrid),      POINTER         :: vgrid_canopy_q_assimi  !< Vertical grid
    TYPE(t_jsb_model),      POINTER         :: model
    TYPE(t_lnd_bgcm_store), POINTER         :: bgcm_store             !< the bgcm store of this tile
    LOGICAL,  DIMENSION(options%nc)         :: l_growing_season       !< growing_season LOGICAL
    LOGICAL,  ALLOCATABLE, DIMENSION(:,:)   :: l_growing_season_cl    !< growing_season LOGICAL, at vgrid_canopy_q_assimi
    INTEGER                                 :: icanopy                !< looping
    INTEGER                                 :: ncanopy                !< number of canopy layers, from vgrid
    INTEGER                                 :: iblk, ic, ics, ice, nc !< dimensions, looping
    REAL(wp)                                :: dtime                  !< timestep
    CHARACTER(len=*), PARAMETER             :: routine = modname//':update_time_average_q_radiation'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_mt2L2D :: veg_pool_mt
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(Q_RAD_)
    dsl4jsb_Def_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: daytime_counter
    dsl4jsb_Real2D_onChunk      :: local_time_day_seconds
    ! Q_PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! 2D Q_RAD_
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc
    dsl4jsb_Real2D_onChunk      :: rfr_ratio_boc_tvegdyn_mavg
    ! 3D Q_RAD_
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_daytime_dacc_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_daytime_dacc_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_daytime_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_daytime_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tcnl_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_shaded_tfrac_mavg_cl
    dsl4jsb_Real3D_onChunk      :: ppfd_sunlit_tfrac_mavg_cl
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_RAD_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)                 RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model                 => Get_model(tile%owner_model_id)
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    ncanopy               =  vgrid_canopy_q_assimi%n_levels
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(Q_RAD_)
    dsl4jsb_Get_memory(Q_PHENO_)

    ALLOCATE(l_growing_season_cl(nc, ncanopy))

    bgcm_store => tile%bgcm_store
    dsl4jsb_Get_mt2L2D(VEG_BGCM_POOL_ID, veg_pool_mt)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)                 ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, daytime_counter)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, local_time_day_seconds)         ! in
    ! Q_PHENO_
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)             ! in
    ! 2D Q_RAD_
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, rfr_ratio_boc)                  ! in
    dsl4jsb_Get_var2D_onChunk(Q_RAD_, rfr_ratio_boc_tvegdyn_mavg)     ! inout
    ! 3D Q_RAD_
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_daytime_dacc_cl)    ! inout
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_daytime_dacc_cl)    ! inout
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_daytime_cl)         ! inout
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_daytime_cl)         ! inout
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_cl)                 ! in
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_cl)                 ! in
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_tcnl_mavg_cl)       ! inout
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_tcnl_mavg_cl)       ! inout
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_shaded_tfrac_mavg_cl)      ! inout
    dsl4jsb_Get_var3D_onChunk(Q_RAD_, ppfd_sunlit_tfrac_mavg_cl)      ! inout
    ! ----------------------------------------------------------------------------------------------------- !

    !>0.9 transform REAL growing_season(:) into LOGICAL l_growing_season(:)/_cl(:,:)
    !>
    DO ic = 1, nc
      IF (growing_season(ic) > test_false_true) THEN
        l_growing_season(ic)      = .TRUE.
        l_growing_season_cl(ic,:) = .TRUE.
      ELSE
        l_growing_season(ic)      = .FALSE.
        l_growing_season_cl(ic,:) = .FALSE.
      END IF
    END DO

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 daytime averages
    !>
    !>

    !>  1.1 accumulate values - if light is available (i.e. daytime)
    !>
    DO ic = 1, nc
      IF (swpar_srf_down(ic) > eps8) THEN
        ! add values
        DO icanopy = 1, ncanopy
          ppfd_sunlit_daytime_dacc_cl(ic,icanopy)  = ppfd_sunlit_daytime_dacc_cl(ic,icanopy) + &
            & ppfd_sunlit_cl(ic,icanopy)
          ppfd_shaded_daytime_dacc_cl(ic,icanopy)  = ppfd_shaded_daytime_dacc_cl(ic,icanopy) + &
            & ppfd_shaded_cl(ic,icanopy)
        END DO
      END IF
    END DO

    !>  1.2 calc daytime averages - at the timestep after local midnight (1st timestep of the new day)
    !>
    !>    calculate the average of the previous day and pass this value to the according variable
    DO ic = 1, nc
      IF (ABS(local_time_day_seconds(ic) - dtime) < eps1) THEN
        ! check if at least one value had been assigned at the current day, avoiding division by zero
        IF (daytime_counter(ic) > eps8) THEN
          ppfd_sunlit_daytime_cl(ic,:) = ppfd_sunlit_daytime_dacc_cl(ic,:) / daytime_counter(ic)
          ppfd_shaded_daytime_cl(ic,:) = ppfd_shaded_daytime_dacc_cl(ic,:) / daytime_counter(ic)
        ELSE
          ppfd_sunlit_daytime_cl(ic,:) = 0.0_wp
          ppfd_shaded_daytime_cl(ic,:) = 0.0_wp
        END IF
        ! zero the accumulation variables after daily average has been calculated
        ppfd_sunlit_daytime_dacc_cl(ic,:) = 0.0_wp
        ppfd_shaded_daytime_dacc_cl(ic,:) = 0.0_wp
      END IF
    END DO

    ! ----------------------------------------------------------------------------------------------------- !
    !>2.0 moving averages
    !>
    !>

    ! docu:
    ! calc_time_mavg(current average, new value, length of avg_period,  !
    !                do_calc=LOGICAL, avg_period_unit='day')            ! OPTIONAL
    !                RETURN(new current average)
    ! the unit of the averaging period is 'day' by default, but can also be 'week' or 'year'

    !>  2.1 tfrac (averages at the time-scale of within-leaf N allocation fractions)
    !>
    ppfd_sunlit_tfrac_mavg_cl(:,:) = calc_time_mavg(ppfd_sunlit_tfrac_mavg_cl(:,:), ppfd_sunlit_daytime_cl(:,:), &
                                                    mavg_period_tfrac, do_calc=l_growing_season_cl(:,:))
    ppfd_shaded_tfrac_mavg_cl(:,:) = calc_time_mavg(ppfd_shaded_tfrac_mavg_cl(:,:), ppfd_shaded_daytime_cl(:,:), &
                                                    mavg_period_tfrac, do_calc=l_growing_season_cl(:,:))

    !>  2.2 tcnl (averages at the time-scale of leaf N allocation fractions)
    !>
    ppfd_sunlit_tcnl_mavg_cl(:,:) = calc_time_mavg(ppfd_sunlit_tcnl_mavg_cl(:,:), ppfd_sunlit_cl(:,:), &
                                                    mavg_period_tcnl, do_calc=l_growing_season_cl(:,:))
    ppfd_shaded_tcnl_mavg_cl(:,:) = calc_time_mavg(ppfd_shaded_tcnl_mavg_cl(:,:), ppfd_shaded_cl(:,:), &
                                                    mavg_period_tcnl, do_calc=l_growing_season_cl(:,:))

    !>  2.3 tvegdyn (averages at the vegetation dynamics time-scale)
    !!
    rfr_ratio_boc_tvegdyn_mavg(:) = calc_time_mavg(rfr_ratio_boc_tvegdyn_mavg(:), rfr_ratio_boc(:), &
      &                               mavg_period_tvegdyn, do_calc=(veg_pool_mt(ix_leaf,ixC,:) > eps8))

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')
  END SUBROUTINE update_time_average_q_radiation

  ! ======================================================================================================= !
  !> aggregate: time moving averages and daytime averages for Q_RAD_
  !>
  SUBROUTINE aggregate_time_average_q_radiation(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_RAD_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_q_radiation'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    model => Get_model(tile%owner_model_id)
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_RAD_)

    dsl4jsb_Aggregate_onChunk(Q_RAD_, rfr_ratio_boc_tvegdyn_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_sunlit_daytime_dacc_cl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_shaded_daytime_dacc_cl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_shaded_daytime_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_sunlit_daytime_cl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_shaded_tcnl_mavg_cl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_sunlit_tcnl_mavg_cl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_shaded_tfrac_mavg_cl     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_sunlit_tfrac_mavg_cl     , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_time_average_q_radiation

#endif
END MODULE mo_q_rad_interface
