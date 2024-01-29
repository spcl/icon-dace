!> vegetation process interface (QUINCY)
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_veg_interface
#ifndef __NO_JSBACH__


  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process
  ! USE mo_veg_turnover,                ONLY: update_veg_turnover
  ! USE mo_veg_dynamics,                ONLY: update_veg_dynamics
  ! USE mo_veg_growth,                  ONLY: update_veg_growth
  ! USE mo_veg_update_pools,            ONLY: update_veg_pools


  ! Use of processes in this module (USE mo_jsb_process_class, ONLY: )
  dsl4jsb_Use_processes VEG_

  ! Use of process configurations (t_PROC_config)
  dsl4jsb_Use_config(VEG_)

  ! Use of process memories (t_PROC_memory)
  dsl4jsb_Use_memory(VEG_)
  

  IMPLICIT NONE
  PRIVATE
!  PUBLIC :: reset_veg_fluxes, update_veg_turnover, update_veg_dynamics, update_veg_growth
!  PUBLIC :: update_veg_pools
  PUBLIC :: Register_veg_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_veg_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_veg_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_reset_veg_fluxes    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_veg_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_veg_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_veg_turnover task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_turnover
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_veg_turnover      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_turnover !< Aggregates computed task variables
  ! END TYPE tsk_update_veg_turnover
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_veg_dynamics task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_dynamics
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_veg_dynamics      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_dynamics !< Aggregates computed task variables
  ! END TYPE tsk_update_veg_dynamics
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_veg_growth task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_growth
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_veg_growth      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_growth !< Aggregates computed task variables
  ! END TYPE tsk_update_veg_growth
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Type definition: update_veg_pools task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! TYPE, EXTENDS(t_jsb_process_task) :: tsk_update_veg_pools
  ! CONTAINS
  !   PROCEDURE, NOPASS :: Integrate => update_veg_pools      !< Advances task computation for one timestep
  !   PROCEDURE, NOPASS :: Aggregate => aggr_update_veg_pools !< Aggregates computed task variables
  ! END TYPE tsk_update_veg_pools  


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_veg_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_veg_fluxes
    PROCEDURE Create_task_reset_veg_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_veg_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_veg_turnover task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_veg_turnover
  !   PROCEDURE Create_task_update_veg_turnover         !< Constructor function for task
  ! END INTERFACE tsk_update_veg_turnover
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_veg_dynamics task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_veg_dynamics
  !   PROCEDURE Create_task_update_veg_dynamics         !< Constructor function for task
  ! END INTERFACE tsk_update_veg_dynamics
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_veg_growth task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_veg_growth
  !   PROCEDURE Create_task_update_veg_growth         !< Constructor function for task
  ! END INTERFACE tsk_update_veg_growth
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor interface: update_veg_pools task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! INTERFACE tsk_update_veg_pools
  !   PROCEDURE Create_task_update_veg_pools         !< Constructor function for task
  ! END INTERFACE tsk_update_veg_pools


  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_veg_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: VEG_ process
  !! 
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_veg_tasks(this, model_id)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id


    CALL this%Register_task(tsk_reset_veg_fluxes(model_id))
    ! CALL this%Register_task(tsk_update_veg_turnover(model_id))
    ! CALL this%Register_task(tsk_update_veg_dynamics(model_id))
    ! CALL this%Register_task(tsk_update_veg_growth(model_id))
    ! CALL this%Register_task(tsk_update_veg_pools(model_id))

  END SUBROUTINE Register_veg_tasks


  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_veg_fluxes task
  !! 
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_veg_fluxes(model_id) RESULT(return_ptr)

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr


    ALLOCATE(tsk_reset_veg_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_veg_fluxes', process_id=VEG_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_veg_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_veg_turnover task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_veg_turnover(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_veg_turnover::return_ptr)
  !   CALL return_ptr%Construct(name='update_veg_turnover', process_id=VEG_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_veg_turnover
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_veg_dynamics task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_veg_dynamics(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_veg_dynamics::return_ptr)
  !   CALL return_ptr%Construct(name='update_veg_dynamics', process_id=VEG_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_veg_dynamics
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_veg_growth task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_veg_growth(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_veg_growth::return_ptr)
  !   CALL return_ptr%Construct(name='update_veg_growth', process_id=VEG_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_veg_growth
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Constructor: update_veg_pools task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! FUNCTION Create_task_update_veg_pools(model_id) RESULT(return_ptr)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   INTEGER,                   INTENT(in) :: model_id
  !   CLASS(t_jsb_process_task), POINTER    :: return_ptr
  ! 
  ! 
  !   ALLOCATE(tsk_update_veg_pools::return_ptr)
  !   CALL return_ptr%Construct(name='update_veg_pools', process_id=VEG_, owner_model_id=model_id)
  ! 
  ! END FUNCTION Create_task_update_veg_pools

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "update": reset_veg_fluxes task
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

  !-----------------------------------------------------------------------------------------------------
  !> Implementation of "aggregate": reset_veg_fluxes task
  !! Actually NOT using this routine, i.e., no var is aggregated here
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE aggregate_reset_veg_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local
    TYPE(t_jsb_model),        POINTER :: model
    CLASS(t_jsb_aggregator),  POINTER :: weighted_by_fract
    INTEGER                           :: iblk, ics, ice, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_veg_fluxes'

    ! Declare process configuration and memory Pointers
    !dsl4jsb_Def_config(VEG_)
    dsl4jsb_Def_memory(VEG_)

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process config
    !dsl4jsb_Get_config(VEG_)
    ! Get process memories
    dsl4jsb_Get_memory(VEG_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! nothing to aggregate, yet

    IF (debug_on() .AND. iblk==1) CALL message(routine, 'Finished.')

  END SUBROUTINE aggregate_reset_veg_fluxes

  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_veg_turnover task
  ! !!  actually NOT using this routine, i.e., no var is aggregated here
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_veg_turnover(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_turnover'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(VEG_)
  !   dsl4jsb_Def_memory(VEG_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(VEG_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(VEG_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  !     ! do not aggregate after this routine, only after update_veg_pools()
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_veg_turnover
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_veg_dynamics task
  ! !!  actually NOT using this routine, i.e., no var is aggregated here
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_veg_dynamics(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_dynamics'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(VEG_)
  !   dsl4jsb_Def_memory(VEG_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(VEG_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(VEG_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  !     ! do not aggregate after this routine, only after update_veg_pools()
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_veg_dynamics
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_veg_growth task
  ! !!  actually NOT using this routine, i.e., no var is aggregated here
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_veg_growth(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_growth'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(VEG_)
  !   dsl4jsb_Def_memory(VEG_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(VEG_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(VEG_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  !     ! do not aggregate after this routine, only after update_veg_pools()
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_veg_growth
  ! 
  ! !-----------------------------------------------------------------------------------------------------
  ! !> Implementation of "aggregate": update_veg_pools task
  ! !! 
  ! !-----------------------------------------------------------------------------------------------------
  ! SUBROUTINE aggr_update_veg_pools(tile, options)
  ! 
  !   IMPLICIT NONE
  !   ! ---------------------------
  !   ! 0.1 InOut
  !   CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
  !   TYPE(t_jsb_task_options),   INTENT(in)    :: options
  !   ! ---------------------------
  !   ! 0.2 Local
  !   TYPE(t_jsb_model),        POINTER         :: model
  !   CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
  !   INTEGER                                   :: iblk, ics, ice, nc
  !   !X if necessary: REAL(wp) :: dtime, steplen
  !   CHARACTER(len=*),         PARAMETER       :: routine = TRIM(modname)//':aggr_update_veg_pools'
  !   ! ---------------------------
  !   ! 0.3 Declare Memory
  !   ! Declare process configuration and memory Pointers
  !   dsl4jsb_Def_config(VEG_)
  !   dsl4jsb_Def_memory(VEG_)
  !   ! Get local variables from options argument
  !   iblk    = options%iblk
  !   ics     = options%ics
  !   ice     = options%ice
  !   nc      = options%nc
  !   !X if necessary: dtime   = options%dtime
  !   !X if necessary: steplen = options%steplen
  !   ! ---------------------------
  !   ! 0.4 Process Activity, Debug Option
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
  !   ! ---------------------------
  !   ! 0.5 Get Memory
  !   ! Get process variables (Set pointers to variables in memory)
  !   model => Get_model(tile%owner_model_id)
  !   ! Get process config
  !   dsl4jsb_Get_config(VEG_)
  !   ! Get process memories
  !   dsl4jsb_Get_memory(VEG_)
  ! 
  ! 
  !   ! ------------------------------------------------------------------------------------------------------------
  !   ! Go aggregate
  !   ! ------------------------------------------------------------------------------------------------------------
  ! 
  ! 
  !   IF (tile%Has_children()) THEN
  !     weighted_by_fract => tile%Get_aggregator("weighted_by_fract")
  ! 
  !     ! *slmdev 
  !     ! the below list is not complete, it is just an initial placeholder for all var that need to be aggregated !
  !     dsl4jsb_Aggregate_onChunk(VEG_, diameter            , weighted_by_fract)
  !     dsl4jsb_Aggregate_onChunk(VEG_, height              , weighted_by_fract)
  !     dsl4jsb_Aggregate_onChunk(VEG_, mean_leaf_age       , weighted_by_fract)
  !     dsl4jsb_Aggregate_onChunk(VEG_, dphi                , weighted_by_fract)
  !     dsl4jsb_Aggregate_onChunk(VEG_, npp                 , weighted_by_fract)
  !   END IF
  ! 
  !   !X Implementation: Start your process scheme here...
  !   !R: Das ist ein Provisorium bisher. Es muss noch festgelegt werden ob und wie die einzelnen Tasks aggregiert werden sollen.
  ! 
  !   IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
  ! 
  ! END SUBROUTINE aggr_update_veg_pools

#endif
END MODULE mo_veg_interface
