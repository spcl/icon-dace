!> Contains the interfaces to the fuel process
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

MODULE mo_fuel_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes FUEL_, CARBON_, PHENO_

  ! Use of process configurations
  dsl4jsb_Use_config(FUEL_)

  ! Use of process memories
  dsl4jsb_Use_memory(FUEL_)
  dsl4jsb_Use_memory(CARBON_)
  dsl4jsb_Use_memory(PHENO_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_fuel_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_interface'

  !> Type definition for fuel task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_fuel
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_fuel    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_fuel !< Aggregates computed task variables
  END TYPE tsk_fuel

  !> Constructor interface for fuel task
  INTERFACE tsk_fuel
    PROCEDURE Create_task_fuel     !< Constructor function for task
  END INTERFACE tsk_fuel

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for fuel task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "fuel"
  !!
  FUNCTION Create_task_fuel(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_fuel::return_ptr)
    CALL return_ptr%Construct(name='fuel', process_id=FUEL_, owner_model_id=model_id)

  END FUNCTION Create_task_fuel

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for fuel process
  !!
  !! @param[in,out] this      Instance of fuel process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_fuel_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_fuel(model_id))

  END SUBROUTINE Register_fuel_tasks


  ! ================================================================================================================================
  !>
  !> Implementation of task fuel
  !! Task "fuel" calculates the fuel available for the fire process
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_fuel(tile, options)

    USE mo_fuel_process,     ONLY: calc_fuel_jsbach
    USE mo_jsb_time,         ONLY: is_newday
    USE mo_jsb_tile_class,   ONLY: t_jsb_tile_abstract

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model

    INTEGER :: iblk, ics, ice

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_fuel'

    ! Declare process memories
    dsl4jsb_Def_memory(CARBON_)
    dsl4jsb_Def_memory(FUEL_)
    dsl4jsb_Def_memory(PHENO_)

    dsl4jsb_Def_config(FUEL_)

    dsl4jsb_Real2D_onChunk :: fuel                   ! Amount of fuel on which fire is based
    !
    dsl4jsb_Real2D_onChunk :: c_acid_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_water_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag1_ta
    dsl4jsb_Real2D_onChunk :: c_acid_ag2_ta
    dsl4jsb_Real2D_onChunk :: c_water_ag2_ta
    dsl4jsb_Real2D_onChunk :: c_ethanol_ag2_ta
    dsl4jsb_Real2D_onChunk :: c_nonsoluble_ag2_ta
    !
    ! Locally allocated vectors
    !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(FUEL_)) RETURN

    IF (.NOT. tile%Has_process_memory(CARBON_)) &
      & CALL finish(TRIM(routine), 'FUEL process depends on CARBON process')

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get model
    model   => Get_model(tile%owner_model_id)

    ! Get process configurations
    dsl4jsb_Get_config(FUEL_)

    ! Get process memories
    dsl4jsb_Get_memory(FUEL_)
    dsl4jsb_Get_memory(CARBON_)
    dsl4jsb_Get_memory(PHENO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag1_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag1_ta)        ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag1_ta)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag1_ta)   ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_acid_ag2_ta)         ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_water_ag2_ta)        ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_ethanol_ag2_ta)      ! in
    dsl4jsb_Get_var2D_onChunk(CARBON_,  c_nonsoluble_ag2_ta)   ! in
    dsl4jsb_Get_var2D_onChunk(FUEL_,    fuel)                  ! out


    !
    ! ---------------------------
    ! Go

    IF (is_newday(options%current_datetime, options%dtime)) THEN ! only once per day

      SELECT CASE (dsl4jsb_Config(FUEL_)%fuel_algorithm)
        CASE (1) !! jsbach algorithm

          CALL calc_fuel_jsbach(       &
            & c_acid_ag1_ta(:),        & ! in
            & c_water_ag1_ta(:),       & ! in
            & c_ethanol_ag1_ta(:),     & ! in
            & c_nonsoluble_ag1_ta(:),  & ! in
            & c_acid_ag2_ta(:),        & ! in
            & c_water_ag2_ta(:),       & ! in
            & c_ethanol_ag2_ta(:),     & ! in
            & c_nonsoluble_ag2_ta(:),  & ! in
            & fuel(:)                  & ! out
            & )

        CASE (2) !! Arora & Boer algorithm
          CALL finish(TRIM(routine),'Arora & Boer fuel algorithm not implemented yet.')
        CASE (3)
          CALL finish(TRIM(routine),'Thornicke/spitfire fuel algorithm not implemented yet.')
        CASE (4) !! Read GFED
        CASE DEFAULT
          CALL finish(TRIM(routine),'Unknown fuel algorithm')
      END SELECT

    ENDIF ! is_newday

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_fuel

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "fuel"
  !!
  !! @param[in,out] tile    Tile for which aggregation of child tiles is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_fuel(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(FUEL_)

    CLASS(t_jsb_aggregator), POINTER          :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_fuel'

    INTEGER  :: iblk, ics, ice

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(FUEL_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(FUEL_, fuel, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_fuel

#endif
END MODULE mo_fuel_interface
