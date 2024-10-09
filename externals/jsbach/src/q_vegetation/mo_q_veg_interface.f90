!> QUINCY vegetation process interface
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
!>#### definition and init of tasks for the vegetation process incl. update and aggregate routines
!>
!> includes plant growth and turnover
!>
MODULE mo_q_veg_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process

  ! "Integrate" routines
  USE mo_q_veg_turnover,              ONLY: update_veg_turnover_real => update_veg_turnover
  USE mo_q_veg_dynamics,              ONLY: update_veg_dynamics_real => update_veg_dynamics
  USE mo_q_veg_growth,                ONLY: update_veg_growth_real => update_veg_growth
  USE mo_q_veg_update_pools,          ONLY: update_veg_pools_real => update_veg_pools
  USE mo_q_veg_products_decay,        ONLY: update_products_decay_real => update_products_decay

  dsl4jsb_Use_processes VEG_, Q_RAD_, SB_
  dsl4jsb_Use_memory(VEG_)
  dsl4jsb_Use_config(VEG_)
  dsl4jsb_Use_memory(Q_RAD_)
  dsl4jsb_Use_memory(SB_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_veg_tasks_quincy

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_veg_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_veg_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_veg_fluxes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_veg_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_turnover task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_turnover
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_turnover     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_turnover  !< Aggregates computed task variables
  END TYPE tsk_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_dynamics task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_dynamics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_dynamics     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_dynamics  !< Aggregates computed task variables
  END TYPE tsk_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_growth task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_growth
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_growth     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_growth  !< Aggregates computed task variables
  END TYPE tsk_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_veg_pools task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_veg_pools
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_veg_pools      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_veg_pools   !< Aggregates computed task variables
  END TYPE tsk_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: time_average_vegetation task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_vegetation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_vegetation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_vegetation !< Aggregates computed task variables
  END TYPE tsk_time_average_vegetation

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_products_decay task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_products_decay
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_products_decay    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_products_decay !< Aggregates computed task variables
  END TYPE tsk_products_decay

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_veg_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_veg_fluxes
    PROCEDURE Create_task_reset_veg_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_veg_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_turnover task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_turnover
    PROCEDURE Create_task_update_veg_turnover         !< Constructor function for task
  END INTERFACE tsk_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_dynamics task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_dynamics
    PROCEDURE Create_task_update_veg_dynamics         !< Constructor function for task
  END INTERFACE tsk_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_growth task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_growth
    PROCEDURE Create_task_update_veg_growth         !< Constructor function for task
  END INTERFACE tsk_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_veg_pools task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_veg_pools
    PROCEDURE Create_task_update_veg_pools         !< Constructor function for task
  END INTERFACE tsk_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_time_average_vegetation task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_time_average_vegetation
    PROCEDURE Create_task_update_time_average_vegetation         !< Constructor function for task
  END INTERFACE tsk_time_average_vegetation

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_products_decay task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_products_decay
    PROCEDURE Create_task_update_products_decay         !< Constructor function for task
  END INTERFACE tsk_products_decay

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_veg_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: VEG_ process
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_veg_tasks_quincy(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_veg_fluxes(model_id))
    CALL this%Register_task(tsk_veg_turnover(model_id))
    CALL this%Register_task(tsk_veg_dynamics(model_id))
    CALL this%Register_task(tsk_veg_growth(model_id))
    CALL this%Register_task(tsk_veg_pools(model_id))
    CALL this%Register_task(tsk_time_average_vegetation(model_id))
    CALL this%Register_task(tsk_products_decay(model_id))

  END SUBROUTINE Register_veg_tasks_quincy

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_veg_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_veg_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_veg_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_veg_fluxes', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_veg_fluxes


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_turnover task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_turnover(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_turnover::return_ptr)
    CALL return_ptr%Construct(name='veg_turnover', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_turnover

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_dynamics task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_dynamics(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_dynamics::return_ptr)
    CALL return_ptr%Construct(name='veg_dynamics', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_dynamics

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_growth task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_growth(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_growth::return_ptr)
    CALL return_ptr%Construct(name='veg_growth', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_growth

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_veg_pools task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_veg_pools(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_veg_pools::return_ptr)
    CALL return_ptr%Construct(name='veg_pools', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: tsk_time_average_vegetation task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_time_average_vegetation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_vegetation::return_ptr)
    CALL return_ptr%Construct(name='tavrg_vegetation', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_time_average_vegetation

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_products_decay task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_products_decay(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_products_decay::return_ptr)
    CALL return_ptr%Construct(name='products_decay', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_update_products_decay

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "update": reset_veg_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE update_reset_veg_fluxes(tile, options)

    USE mo_veg_util, ONLY: reset_veg_fluxes

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_veg_fluxes'

    iblk    = options%iblk

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL reset_veg_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_reset_veg_fluxes

  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrappers for update routines from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_veg_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_pools_real(tile, options)

  END SUBROUTINE update_veg_pools

  SUBROUTINE update_veg_growth(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_growth_real(tile, options)

  END SUBROUTINE update_veg_growth

  SUBROUTINE update_veg_dynamics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_dynamics_real(tile, options)

  END SUBROUTINE update_veg_dynamics

  SUBROUTINE update_veg_turnover(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_veg_turnover_real(tile, options)

  END SUBROUTINE update_veg_turnover

  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrappers for update routines from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_products_decay(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_products_decay_real(tile, options)

  END SUBROUTINE update_products_decay

  ! ======================================================================================================= !
  !>Implementation of "aggregate": reset_veg_fluxes task
  !>
  SUBROUTINE aggregate_reset_veg_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_veg_fluxes'
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
    dsl4jsb_Get_memory(VEG_)

    ! VEG_ 2D
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_pot             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_c13             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_c14             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_respiration                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_respiration_c13            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_respiration_c14            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_transform_respiration           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation_respiration            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_processing_respiration          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_processing_respiration_c13      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_processing_respiration_c14      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp                               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_c13                           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_c14                           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_growth                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_nh4                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_nh4_n15                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_no3                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_no3_n15                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation_n15                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_po4                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_leaf_n                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_leaf_n15                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_leaf_p                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_fine_root_n             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_fine_root_n15           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_fine_root_p             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_heart_wood_n            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_heart_wood_n15          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_heart_wood_p            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, delta_dens_ind                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_transpiration                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_biosphere_production          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, biological_n_fixation             , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_reset_veg_fluxes


  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_turnover task
  !>
  SUBROUTINE aggregate_veg_turnover(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_turnover'
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
    dsl4jsb_Get_memory(VEG_)

    ! VEG_ 2D
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_leaf_n                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_fine_root_n           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_leaf_p                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_fine_root_p           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_leaf_n15              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_fine_root_n15         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_heart_wood_n          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_heart_wood_n15        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, recycling_heart_wood_p          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_growth                      , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_veg_turnover

  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_dynamics task
  !>
  SUBROUTINE aggregate_veg_dynamics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_dynamics'
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
    dsl4jsb_Get_memory(VEG_)

    ! VEG_ 2D
    dsl4jsb_Aggregate_onChunk(VEG_, do_cohort_harvest           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, cohort_age                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, mortality_rate              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, delta_dens_ind              , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_veg_dynamics

  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_growth task
  !>
  SUBROUTINE aggregate_veg_growth(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_growth'
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
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_memory(SB_)

    ! VEG_ 2D
    dsl4jsb_Aggregate_onChunk(VEG_, target_cn_leaf                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, leaf_cn_direction                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_npp                          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_transpiration                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, kstar_labile                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fmaint_rate_root                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_n_act                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_p_act                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, cost_n_uptake_root                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_np_leaf                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_cn_fine_root               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_cn_coarse_root             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_cn_sap_wood                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_cn_fruit                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_np_fine_root               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_np_coarse_root             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_np_sap_wood                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_np_fruit                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, leaf2sapwood_mass_ratio           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, leaf2root_mass_ratio              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, target_lai                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, k1_opt                            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, k2_opt                            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_p                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_n                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_nh4                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_nh4_n15                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_no3                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_no3_n15                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_po4                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, f_n_demand                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, f_p_demand                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_transform_respiration           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation_respiration            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_processing_respiration          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_c13             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_c14             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_respiration_c13            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_respiration_c14            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_processing_respiration_c13      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_processing_respiration_c14      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_pot             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_respiration                , weighted_by_fract)
    ! SB_ 3D
    dsl4jsb_Aggregate_onChunk(SB_, plant_uptake_nh4_sl                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, plant_uptake_nh4_n15_sl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, plant_uptake_no3_sl                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, plant_uptake_no3_n15_sl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, plant_uptake_po4_sl                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tlabile_mavg_sl       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_p_tlabile_mavg_sl       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_c_tmyc_mavg_sl          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tmyc_mavg_sl          , weighted_by_fract)


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_veg_growth

  ! ======================================================================================================= !
  !>Implementation of "aggregate": update_veg_pools task
  !>
  SUBROUTINE aggregate_veg_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(Q_RAD_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_veg_pools'
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
    dsl4jsb_Get_memory(VEG_)
    dsl4jsb_Get_config(VEG_)
    dsl4jsb_Get_memory(Q_RAD_)

    ! VEG_ 2D
    dsl4jsb_Aggregate_onChunk(VEG_, diameter                            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, height                              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, dens_ind                            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, lai                                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_growth                          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, mean_leaf_age                       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fract_fpc                           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, blended_height                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, dphi                                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp                                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_c13                             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_c14                             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation_n15                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_biosphere_production            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, biological_n_fixation               , weighted_by_fract)
    ! VEG_ bgcm pools
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_c            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_n            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_p            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_c13          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_c14          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_total_n15          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_leaf_c             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_leaf_n             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_wood_c             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_wood_n             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_fine_root_c        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_pool_fine_root_n        , weighted_by_fract)

    IF(dsl4jsb_Config(VEG_)%l_use_product_pools) THEN
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_c        , weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(VEG_, veg_products_total_n        , weighted_by_fract)
    ENDIF
    ! VEG_ bgcm fluxes
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_c          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_n          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_p          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_c13        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_c14        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_growth_total_n15        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_c      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_n      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_p      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_c13    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_c14    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, veg_litterfall_total_n15    , weighted_by_fract)
    ! VEG_ 3D
    dsl4jsb_Aggregate_onChunk(VEG_, root_fraction_sl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, leaf_nitrogen_cl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fn_rub_cl                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fn_et_cl                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fn_pepc_cl                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fn_chl_cl                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fn_oth_cl                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, lai_cl                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, cumm_lai_cl                 , weighted_by_fract)
    ! Q_RAD_
    dsl4jsb_Aggregate_onChunk(Q_RAD_, rfr_ratio_boc_tvegdyn_mavg  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_sunlit_tfrac_mavg_cl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_RAD_, ppfd_shaded_tfrac_mavg_cl   , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_veg_pools


  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for VEG_
  !>
  SUBROUTINE update_time_average_vegetation(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_veg_util,            ONLY: calculate_time_average_vegetation
    dsl4jsb_Use_processes VEG_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_vegetation'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(VEG_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)           RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL calculate_time_average_vegetation(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_vegetation


  ! ======================================================================================================= !
  !> Implementation of "aggregate": time_average_vegetation task
  !>
  SUBROUTINE aggregate_time_average_vegetation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_vegetation'
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
    dsl4jsb_Get_memory(VEG_)

    ! VEG_ 2D
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_daytime_dacc                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_daytime_dacc            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_daytime_dacc     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_daytime_dacc                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps_daytime_dacc      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_jmax_opt_daytime_dacc           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_daytime                     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_daytime                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_daytime          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_daytime                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps_daytime           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_jmax_opt_daytime                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, gpp_tlabile_mavg                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, maint_respiration_tlabile_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_n_tlabile_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_p_tlabile_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tphen_mavg                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_tuptake_mavg                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, demand_uptake_n_tuptake_mavg      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, demand_uptake_p_tuptake_mavg      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_n_tuptake_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_req_p_tuptake_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tfrac_mavg                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_tfrac_mavg              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_tfrac_mavg       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_tfrac_mavg                     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, beta_sinklim_ps_tfrac_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tcnl_mavg                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, press_srf_tcnl_mavg               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, co2_mixing_ratio_tcnl_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, ga_tcnl_mavg                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, uptake_n_tcnl_mavg                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fmaint_rate_tcnl_mavg             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_cn_tcnl_mavg               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_np_tcnl_mavg               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_tcnl_mavg                     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_carbon_tcnl_mavg           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_nitrogen_tcnl_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_phosphorus_tcnl_mavg       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, npp_talloc_mavg                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, n_fixation_talloc_mavg            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_npp_talloc_mavg              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, transpiration_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_transpiration_talloc_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, dphi_talloc_mavg                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_n_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_p_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_carbon_talloc_mavg         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_nitrogen_talloc_mavg       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, labile_phosphorus_talloc_mavg     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, reserve_carbon_talloc_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, reserve_nitrogen_talloc_mavg      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, reserve_phosphorus_talloc_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, w_root_lim_talloc_mavg            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_cn_talloc_mavg             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_cp_talloc_mavg             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, growth_np_talloc_mavg             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, exudation_c_tmyc_mavg             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, leaf2root_troot_mavg              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_npp_troot_mavg               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_n_pot_troot_mavg      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_p_pot_troot_mavg      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fmaint_rate_troot_mavg            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, an_boc_tvegdyn_mavg               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, net_growth_tvegdyn_mavg           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, lai_tvegdyn_mavg                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_tacclim_mavg                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_soil_root_tacclim_mavg          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_week_mavg                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_air_month_mavg                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, t_jmax_opt_mavg                   , weighted_by_fract)
    ! VEG_ 3D
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_daytime_dacc_cl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_daytime_cl           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_tcnl_mavg_cl         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, fleaf_sunlit_tfrac_mavg_cl        , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_vegetation

  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for products decay task
  !
  SUBROUTINE aggregate_products_decay(tile, options)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_products_decay'
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
    dsl4jsb_Get_memory(VEG_)

    ! VEG_ 2D
    !TODO: implement -- only required with IQ
    !dsl4jsb_Aggregate_onChunk(VEG_, XXX , weighted_by_fract)

  END SUBROUTINE aggregate_products_decay

#endif
END MODULE mo_q_veg_interface
