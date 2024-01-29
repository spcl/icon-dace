!> interface to the tcq process (test conserved quantities)
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
!>#### Contains the interfaces to the tcq process which is a small test process containing a couple of conserved quantities 
!>
!> Note: 
!> 
!>       1. This module takes advantage of the handling of the integrate and aggregate routines, which both
!>          call the same routines in this module (integrate = aggregate). 
!>          Thereby each tile is processed by the same routine for the tasks of this process.
!>
!>       2. currently only 2D variables are dealt with -- if also 3D variables need to be conserved,
!>          these will anyway have to specially be dealt with in the concept and thus need to be 
!>          identified with an own CQT e.g. distinguish WATER_2D_CQ_TYPE and WATER_3D_CQ_TYPE?
!>          Also, integration of pool-structure might need some additional work.
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_tcq_interface
#ifndef __NO_JSBACH__

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model    
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes TCQ_, PHENO_

  ! Use of process configurations
  dsl4jsb_Use_config(TCQ_)

  ! Use of process memories
  dsl4jsb_Use_memory(TCQ_)
  dsl4jsb_Use_memory(PHENO_)

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Register_tcq_tasks, tcq_transfer_from_active_to_passive_vars_onChunk

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tcq_interface'
  
  !> Type definition for tcq task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_tcq
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_tcq     !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_tcq  !< Aggregates computed task variables
  END TYPE tsk_tcq
  
  !> Constructor interface for tcq task
  INTERFACE tsk_tcq
    PROCEDURE Create_task_tcq                        !< Constructor function for task
  END INTERFACE tsk_tcq
  
CONTAINS

  ! ====================================================================================================== !
  !
  !> Constructor for tcq task
  !
  FUNCTION Create_task_tcq(model_id) RESULT(return_ptr)

    ! -------------------------------------------------------------------------------------------------- !
    INTEGER,                   INTENT(in) :: model_id    !< Model id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr  !< Instance of process task "tcq"
    ! -------------------------------------------------------------------------------------------------- !

    ALLOCATE(tsk_tcq::return_ptr)
    CALL return_ptr%Construct(name='tcq', process_id=TCQ_, owner_model_id=model_id)

  END FUNCTION Create_task_tcq
    
  ! ====================================================================================================== !
  !
  !> Register tasks for tcq process
  !
  SUBROUTINE Register_tcq_tasks(this, model_id)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_process), INTENT(inout) :: this        !< Instance of tcq process class
    INTEGER,               INTENT(in)   :: model_id    !< Model id
    ! -------------------------------------------------------------------------------------------------- !

    CALL this%Register_task(tsk_tcq(model_id))
    
  END SUBROUTINE Register_tcq_tasks

  ! ====================================================================================================== !
  !
  !> Implementation of "update" for task "tcq"
  !>
  !> Calculates a ta variable for each tcq variable
  !
  SUBROUTINE update_tcq(tile, options)

    IMPLICIT NONE

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER          :: model

    INTEGER  :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_tcq'

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_memory(TCQ_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real2D_onChunk ::  a_dead_c_tcq
    dsl4jsb_Real2D_onChunk ::  another_dead_c_tcq
    dsl4jsb_Real2D_onChunk ::  a_veg_c_tcq
    dsl4jsb_Real2D_onChunk ::  an_implicit_scaling_tcq

    dsl4jsb_Real2D_onChunk ::  a_dead_c_ta_tcq
    dsl4jsb_Real2D_onChunk ::  another_dead_c_ta_tcq
    dsl4jsb_Real2D_onChunk ::  a_veg_c_ta_tcq
    dsl4jsb_Real2D_onChunk ::  an_implicit_scaling_ta_tcq

    dsl4jsb_Real2D_onChunk ::  veg_fract_correction
    dsl4jsb_Real2D_onChunk ::  fract_fpc_max    
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    
    model => Get_model(tile%owner_model_id)

    ! If process is not active on this tile, do nothing
    IF (.NOT. tile%Is_process_active(TCQ_)) RETURN

    dsl4jsb_Get_memory(TCQ_)
    dsl4jsb_Get_memory(PHENO_)

    dsl4jsb_Get_var2D_onChunk(TCQ_,  a_dead_c_tcq )
    dsl4jsb_Get_var2D_onChunk(TCQ_,  another_dead_c_tcq )
    dsl4jsb_Get_var2D_onChunk(TCQ_,  a_veg_c_tcq )
    dsl4jsb_Get_var2D_onChunk(TCQ_,  an_implicit_scaling_tcq )

    dsl4jsb_Get_var2D_onChunk(TCQ_,  a_dead_c_ta_tcq )
    dsl4jsb_Get_var2D_onChunk(TCQ_,  another_dead_c_ta_tcq )
    dsl4jsb_Get_var2D_onChunk(TCQ_,  a_veg_c_ta_tcq )
    dsl4jsb_Get_var2D_onChunk(TCQ_,  an_implicit_scaling_ta_tcq )

    dsl4jsb_Get_var2D_onChunk(PHENO_,  fract_fpc_max )
    dsl4jsb_Get_var2D_onChunk(PHENO_,  veg_fract_correction )

    a_veg_c_ta_tcq = a_veg_c_tcq * veg_fract_correction * fract_fpc_max
    a_dead_c_ta_tcq = a_dead_c_tcq * veg_fract_correction * fract_fpc_max
    another_dead_c_ta_tcq = another_dead_c_tcq * veg_fract_correction * fract_fpc_max
    an_implicit_scaling_ta_tcq = an_implicit_scaling_tcq * veg_fract_correction * fract_fpc_max
    
  END SUBROUTINE update_tcq
  
  ! ====================================================================================================== !
  !
  !> Implementation of "aggregate" for task "tcq"
  !>
  !> Aggregates ta variables associated with the tcq variables
  !
  SUBROUTINE aggregate_tcq(tile, options)

    USE mo_jsb_time,          ONLY: is_newday

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile     !< Tile for which routine is executed
    TYPE(t_jsb_task_options),   INTENT(in)    :: options  !< Additional run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    !dsl4jsb_Def_config(TCQ_)
    dsl4jsb_Def_memory(TCQ_)

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice
    REAL(wp) :: new_sum(options%nc)
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_tcq'
    
    dsl4jsb_Real2D_onChunk ::  old_sum
    dsl4jsb_Real2D_onChunk ::  conservation_test_field
    dsl4jsb_Real2D_onChunk ::  a_dead_c_ta_tcq
    dsl4jsb_Real2D_onChunk ::  another_dead_c_ta_tcq
    dsl4jsb_Real2D_onChunk ::  a_veg_c_ta_tcq
    dsl4jsb_Real2D_onChunk ::  an_implicit_scaling_ta_tcq
    ! -------------------------------------------------------------------------------------------------- !

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    
    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    
    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(TCQ_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(TCQ_, a_dead_c_ta_tcq, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TCQ_, a_veg_c_ta_tcq, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TCQ_, another_dead_c_ta_tcq, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(TCQ_, an_implicit_scaling_ta_tcq, weighted_by_fract)

    IF(is_newday(options%current_datetime,options%dtime)) THEN
      IF ( tile%name .EQ. 'veg') THEN
        dsl4jsb_Get_var2D_onChunk(TCQ_,  old_sum)
        dsl4jsb_Get_var2D_onChunk(TCQ_,  conservation_test_field)
        dsl4jsb_Get_var2D_onChunk(TCQ_,  a_veg_c_ta_tcq)
        dsl4jsb_Get_var2D_onChunk(TCQ_,  a_dead_c_ta_tcq)
        dsl4jsb_Get_var2D_onChunk(TCQ_,  another_dead_c_ta_tcq)
        dsl4jsb_Get_var2D_onChunk(TCQ_,  an_implicit_scaling_ta_tcq)
        new_sum = a_veg_c_ta_tcq + a_dead_c_ta_tcq + another_dead_c_ta_tcq + an_implicit_scaling_ta_tcq
        
        conservation_test_field = new_sum - old_sum
        old_sum = new_sum
      ENDIF
    ENDIF


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_tcq


  ! ====================================================================================================== !
  !
  !> Relocates content of own active to own passive vars that are part of lcc_relocations
  !
  SUBROUTINE tcq_transfer_from_active_to_passive_vars_onChunk(lcc_relocations, tile, i_tile, options)
    USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc
    USE mo_jsb_task_class,      ONLY: t_jsb_task_options
    USE mo_jsb_lcc_class,       ONLY: transfer_from_active_to_passive_var_onChunk
    USE mo_tcq_constants,       ONLY: tcq_potential_active_vars, tcq_required_passive_vars

    IMPLICIT NONE
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lcc_proc),       INTENT(INOUT) :: lcc_relocations
        !< lcc structure of the calling lcc process for distributing matter from an active to a passive var
    CLASS(t_jsb_tile_abstract), INTENT(in)    :: tile !< tile for which the relocation is be conducted
    INTEGER,                    INTENT(in)    :: i_tile !< index of the tile in lcc structure
    TYPE(t_jsb_task_options),   INTENT(IN)    :: options !< run-time parameters
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':tcq_transfer_from_active_to_passive_vars_onChunk'

    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on() .AND. options%iblk == 1) CALL message(TRIM(routine), 'For '//TRIM(lcc_relocations%name)//' ...')

    ! for tcq the transfer is simple - just put from a_veg to a_dead...
    CALL transfer_from_active_to_passive_var_onChunk(                            &
      & lcc_relocations, i_tile, options, tcq_potential_active_vars(1), tcq_required_passive_vars(1))

  END SUBROUTINE tcq_transfer_from_active_to_passive_vars_onChunk


#endif
END MODULE mo_tcq_interface
