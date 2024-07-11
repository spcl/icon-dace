!> QUINCY soil-biogeochemistry process interface
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
!>#### definition and init of tasks for the soil-biogeochemistry process incl. update and aggregate routines
!>
MODULE mo_q_sb_interface
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
  USE mo_q_sb_rate_modifier,          ONLY: update_sb_rate_modifier_real => update_sb_rate_modifier
  USE mo_q_sb_ssm_main,               ONLY: update_sb_simple_model_real => update_sb_simple_model
  USE mo_q_sb_jsm_main,               ONLY: update_sb_jsm_real => update_sb_jsm
  USE mo_q_sb_update_pools,           ONLY: update_sb_pools_real => update_sb_pools

  dsl4jsb_Use_processes SB_, VEG_
  dsl4jsb_Use_config(SB_)
  dsl4jsb_Use_memory(SB_)
  dsl4jsb_Use_memory(VEG_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_sb_tasks_quincy

  ! ======================================================================================================= !
  !>Type definition: reset_sb_fluxes task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_sb_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_sb_fluxes  !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_sb_fluxes    !< Aggregates computed task variables
  END TYPE tsk_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Type definition: sb_rate_modifier task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_rate_modifier
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_rate_modifier       !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_rate_modifier    !< Aggregates computed task variables
  END TYPE tsk_sb_rate_modifier

  ! ======================================================================================================= !
  !>Type definition: sb_simple_model task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_simple_model
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_simple_model        !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_simple_model     !< Aggregates computed task variables
  END TYPE tsk_sb_simple_model

  ! ======================================================================================================= !
  !>Type definition: sb_jsm task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_jsm
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_jsm      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_jsm   !< Aggregates computed task variables
  END TYPE tsk_sb_jsm

  ! ======================================================================================================= !
  !>Type definition: sb_pools task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_sb_pools
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_sb_pools         !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_sb_pools      !< Aggregates computed task variables
  END TYPE tsk_sb_pools

  ! ======================================================================================================= !
  !>Type definition: time_average_soilbiogeochemistry task
  !>
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_time_average_soilbiogeochemistry
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_time_average_soilbiogeochemistry    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_time_average_soilbiogeochemistry !< Aggregates computed task variables
  END TYPE tsk_time_average_soilbiogeochemistry

  ! ======================================================================================================= !
  !>Constructor interface: reset_sb_fluxes task
  !>
  INTERFACE tsk_reset_sb_fluxes
    PROCEDURE Create_task_reset_sb_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Constructor interface: sb_rate_modifier task
  !>
  INTERFACE tsk_sb_rate_modifier
    PROCEDURE Create_task_sb_rate_modifier         !< Constructor function for task
  END INTERFACE tsk_sb_rate_modifier

  ! ======================================================================================================= !
  !>Constructor interface: sb_simple_model task
  !>
  INTERFACE tsk_sb_simple_model
    PROCEDURE Create_task_sb_simple_model         !< Constructor function for task
  END INTERFACE tsk_sb_simple_model

  ! ======================================================================================================= !
  !>Constructor interface: update_sb_jsm task
  !>
  INTERFACE tsk_sb_jsm
    PROCEDURE Create_task_sb_jsm         !< Constructor function for task
  END INTERFACE tsk_sb_jsm

  ! ======================================================================================================= !
  !>Constructor interface: sb_pools task
  !>
  INTERFACE tsk_sb_pools
    PROCEDURE Create_task_sb_pools         !< Constructor function for task
  END INTERFACE tsk_sb_pools

  ! ======================================================================================================= !
  !>Constructor interface: time_average_soilbiogeochemistry task
  !>
  INTERFACE tsk_time_average_soilbiogeochemistry
    PROCEDURE Create_task_update_time_average_soilbiogeochemistry         !< Constructor function for task
  END INTERFACE tsk_time_average_soilbiogeochemistry

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_q_sb_interface'

CONTAINS

  ! ======================================================================================================= !
  !>Register tasks: SB_ process
  !>
  SUBROUTINE Register_sb_tasks_quincy(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_sb_fluxes(model_id))
    CALL this%Register_task(tsk_sb_rate_modifier(model_id))
    CALL this%Register_task(tsk_sb_simple_model(model_id))
    CALL this%Register_task(tsk_sb_jsm(model_id))
    CALL this%Register_task(tsk_sb_pools(model_id))
    CALL this%Register_task(tsk_time_average_soilbiogeochemistry(model_id))
  END SUBROUTINE Register_sb_tasks_quincy

  ! ======================================================================================================= !
  !>Constructor: reset_sb_fluxes task
  !>
  FUNCTION Create_task_reset_sb_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_sb_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_sb_fluxes', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Constructor: sb_rate_modifier task
  !>
  FUNCTION Create_task_sb_rate_modifier(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_rate_modifier::return_ptr)
    CALL return_ptr%Construct(name='sb_rate_modifier', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_rate_modifier

  ! ======================================================================================================= !
  !>Constructor: sb_simple_model task
  !>
  FUNCTION Create_task_sb_simple_model(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_simple_model::return_ptr)
    CALL return_ptr%Construct(name='sb_simple_model', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_simple_model

  ! ======================================================================================================= !
  !>Constructor: sb_jsm task
  !>
  FUNCTION Create_task_sb_jsm(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_jsm::return_ptr)
    CALL return_ptr%Construct(name='sb_jsm', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_jsm

  ! ======================================================================================================= !
  !>Constructor: sb_pools task
  !>
  FUNCTION Create_task_sb_pools(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_sb_pools::return_ptr)
    CALL return_ptr%Construct(name='sb_pools', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_sb_pools

  ! ======================================================================================================= !
  !>Constructor: tsk_time_average_soilbiogeochemistry task
  !>
  FUNCTION Create_task_update_time_average_soilbiogeochemistry(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_time_average_soilbiogeochemistry::return_ptr)
    CALL return_ptr%Construct(name='tavrg_soilbiogeochemistry', process_id=SB_, owner_model_id=model_id)

  END FUNCTION Create_task_update_time_average_soilbiogeochemistry

  ! ------------------------------------------------------------------------------------------------------- !
  ! update_* routines
  ! ------------------------------------------------------------------------------------------------------- !

  ! ======================================================================================================= !
  !>Implementation of "integrate": reset_sb_fluxes task
  !>
  SUBROUTINE update_reset_sb_fluxes(tile, options)

    USE mo_sb_util, ONLY: real_reset_sb_fluxes => reset_sb_fluxes
    dsl4jsb_Use_processes SB_

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    INTEGER :: iblk
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_reset_sb_fluxes'

    iblk = options%iblk
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL real_reset_sb_fluxes(tile, options)

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_reset_sb_fluxes

  ! ======================================================================================================= !
  !>calculate time moving averages and daytime averages for SB_
  !>
  SUBROUTINE update_time_average_soilbiogeochemistry(tile, options)

    USE mo_jsb_control,         ONLY: debug_on
    USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_sb_util,             ONLY: calculate_time_average_soilbiogeochemistry
    dsl4jsb_Use_processes SB_
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_time_average_soilbiogeochemistry'
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(SB_)) RETURN
    IF (tile%lcts(1)%lib_id == 0)          RETURN  ! only if the present tile is a pft
    IF (debug_on() .AND. options%iblk == 1) &
      & CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    CALL calculate_time_average_soilbiogeochemistry(tile, options)

    IF (debug_on() .AND. options%iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE update_time_average_soilbiogeochemistry

  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrappers for update routines from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_sb_rate_modifier(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_rate_modifier_real(tile, options)

  END SUBROUTINE update_sb_rate_modifier

  SUBROUTINE update_sb_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_pools_real(tile, options)

  END SUBROUTINE update_sb_pools

  SUBROUTINE update_sb_jsm(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_jsm_real(tile, options)

  END SUBROUTINE update_sb_jsm

  SUBROUTINE update_sb_simple_model(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_sb_simple_model_real(tile, options)

  END SUBROUTINE update_sb_simple_model

  ! ------------------------------------------------------------------------------------------------------- !
  ! aggregate_* routines
  ! ------------------------------------------------------------------------------------------------------- !

  ! ======================================================================================================= !
  !>Implementation of "aggregate": reset_sb_fluxes task
  !>
  SUBROUTINE aggregate_reset_sb_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_sb_fluxes'
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
    dsl4jsb_Get_memory(SB_)

    ! nothing to aggregate, yet

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE aggregate_reset_sb_fluxes

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_rate_modifier task
  !>
  SUBROUTINE aggregate_sb_rate_modifier(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_rate_modifier'
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
    dsl4jsb_Get_memory(SB_)

    dsl4jsb_Aggregate_onChunk(SB_, rtm_decomposition              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_depolymerisation           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_mic_uptake                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_plant_uptake               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_sorption                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_desorption                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_hsc                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_nitrification              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_denitrification            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rtm_gasdiffusion               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_decomposition              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_depolymerisation           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_mic_uptake                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_plant_uptake               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_sorption                   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_desorption                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_hsc                        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_nitrification              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, rmm_gasdiffusion               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, anaerobic_volume_fraction_sl   , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_rate_modifier

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_simple_model task
  !>
  SUBROUTINE aggregate_sb_simple_model(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_simple_model'
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
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)

    ! VEG 2D
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_n_pot               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_, unit_uptake_p_pot               , weighted_by_fract)
    ! SB_ 2D
    dsl4jsb_Aggregate_onChunk(SB_,  emission_noy                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  emission_n2o                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  emission_n2                     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  emission_noy_n15                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  emission_n2o_n15                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  emission_n2_n15                 , weighted_by_fract)
    ! SB_ 3D
    dsl4jsb_Aggregate_onChunk(SB_,  het_respiration                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  het_respiration_c13             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  het_respiration_c14             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  plant_uptake_nh4_sl             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  plant_uptake_no3_sl             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  plant_uptake_po4_sl             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  asymb_n_fixation                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  asymb_n_fixation_n15            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  k_bioturb                       , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_simple_model

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_jsm task
  !>
  SUBROUTINE aggregate_sb_jsm(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_jsm'
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
    dsl4jsb_Get_memory(SB_)

    ! @TODO  add missing variables for aggregation

    ! SB_ 3D
    dsl4jsb_Aggregate_onChunk(SB_,  asymb_n_fixation                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  asymb_n_fixation_n15            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_,  k_bioturb                       , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_jsm

  ! ======================================================================================================= !
  !>Implementation of "aggregate": sb_pools task
  !>
  SUBROUTINE aggregate_sb_pools(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_sb_pools'
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
    dsl4jsb_Get_memory(SB_)
    dsl4jsb_Get_memory(VEG_)

    ! fluxes
    dsl4jsb_Aggregate_onChunk(SB_, het_respiration              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, het_respiration_c13          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, het_respiration_c14          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_respiration              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_respiration_c13          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_respiration_c14          , weighted_by_fract)
    ! 2D
    dsl4jsb_Aggregate_onChunk(SB_, ecosystem_total_n_loss       , weighted_by_fract)
    ! pools
    dsl4jsb_Aggregate_onChunk(SB_, total_flux_carbon            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, total_flux_nitrogen          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, total_flux_phosphorus        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, total_flux_carbon13          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, total_flux_carbon14          , weighted_by_fract)
    ! 3D
    dsl4jsb_Aggregate_onChunk(SB_, noy                          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, noy_n15                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, n2o                          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, n2o_n15                      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, n2                           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, n2_n15                       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, total_soil_n                 , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, total_soil_inorg_n           , weighted_by_fract)
    ! pools
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_c              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_n              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_p              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_c13            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_c14            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_n15            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_n          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_p          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c13        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c14        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_n15        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_n       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_p       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c13     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c14     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_n15     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, sb_pool_woody_litter_c       , weighted_by_fract)

    ! VEG_ fluxes
    dsl4jsb_Aggregate_onChunk(VEG_,   net_biosphere_production_l2aveghlp  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_,   biological_n_fixation_l2aveghlp     , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_sb_pools

  ! ======================================================================================================= !
  !>Implementation of "aggregate": time_average_assimilation task
  !>
  SUBROUTINE aggregate_time_average_soilbiogeochemistry(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SB_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_time_average_soilbiogeochemistry'
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
    dsl4jsb_Get_memory(SB_)

    dsl4jsb_Aggregate_onChunk(SB_, microbial_cue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, microbial_nue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, microbial_pue_eff_tmic_mavg    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_c_tmyc_mavg_sl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tmyc_mavg_sl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_n_tlabile_mavg_sl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, myc_export_p_tlabile_mavg_sl   , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, dom_cn_mavg                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, dom_cp_mavg                    , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_c_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_n_mavg        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SB_, enzyme_frac_poly_p_mavg        , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c            , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_n            , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_p            , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c13          , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_c14          , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_som_n15          , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c         , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_n         , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_p         , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c13       , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_c14       , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(SB_, sb_pool_total_litter_n15       , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_time_average_soilbiogeochemistry

#endif
END MODULE mo_q_sb_interface
