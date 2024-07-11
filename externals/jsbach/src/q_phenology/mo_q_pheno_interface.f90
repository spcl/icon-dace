!> QUINCY phenology process interface
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
!>#### definition and init of tasks for the phenology process incl. update and aggregate routines
!>
MODULE mo_q_pheno_interface
#ifndef __NO_QUINCY__

  USE mo_jsb_control,        ONLY: debug_on
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, finish
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_q_pheno_process,    ONLY: update_q_phenology_real => update_q_phenology

  dsl4jsb_Use_processes Q_PHENO_, VEG_
  dsl4jsb_Use_config(Q_PHENO_)
  dsl4jsb_Use_memory(Q_PHENO_)
  dsl4jsb_Use_memory(VEG_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_q_pheno_tasks

  ! ======================================================================================================= !
  !> Type definition for the phenology task
  !>
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_q_phenology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_q_phenology    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_q_phenology !< Aggregates computed task variables
  END TYPE   tsk_q_phenology
  !> Constructor interface for update_q_phenology task
  !>
  INTERFACE tsk_q_phenology
    PROCEDURE Create_task_q_phenology         !< Constructor function for task
  END INTERFACE tsk_q_phenology

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_pheno_interface'

CONTAINS

  ! ======================================================================================================= !
  !> Constructor for phenology task
  !>
  FUNCTION Create_task_q_phenology(model_id) RESULT(return_ptr)
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_q_phenology::return_ptr)
    CALL return_ptr%Construct(name='q_phenology', process_id=Q_PHENO_, owner_model_id=model_id)
  END FUNCTION Create_task_q_phenology

  ! ======================================================================================================= !
  !> Register tasks for pheno process
  !>
  SUBROUTINE Register_q_pheno_tasks(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_q_phenology(model_id))
  END SUBROUTINE Register_q_pheno_tasks

  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrapper for update routine from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_q_phenology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_q_phenology_real(tile, options)

  END SUBROUTINE update_q_phenology

  ! ======================================================================================================= !
  !> Implementation of "aggregate" of the "phenology" task
  !>
  SUBROUTINE aggregate_q_phenology(tile, options)
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(Q_PHENO_)
    dsl4jsb_Def_memory(VEG_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_q_phenology'
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    ! ----------------------------------------------------------------------------------------------------- !
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(Q_PHENO_)
    dsl4jsb_Get_memory(VEG_)

    ! Q_PHENO_
    dsl4jsb_Aggregate_onChunk(Q_PHENO_, growing_season              ,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_PHENO_, gdd                         ,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(Q_PHENO_, nd_dormance                 ,    weighted_by_fract)
    ! VEG_
    dsl4jsb_Aggregate_onChunk(VEG_,   growth_req_n_tlabile_mavg   ,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(VEG_,   growth_req_p_tlabile_mavg   ,    weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  END SUBROUTINE aggregate_q_phenology

#endif
END MODULE mo_q_pheno_interface
