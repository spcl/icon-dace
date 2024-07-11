!> QUINCY soil-physics process interface
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
!>#### definition and init of tasks for the soil-physics-quincy process incl. update and aggregate routines
!>
MODULE mo_q_spq_interface
#ifndef __NO_QUINCY__

  USE mo_kind,                        ONLY: wp
  USE mo_jsb_control,                 ONLY: debug_on
  USE mo_exception,                   ONLY: message, finish
  USE mo_jsb_class,                   ONLY: Get_model
  USE mo_jsb_model_class,             ONLY: t_jsb_model
  USE mo_jsb_tile_class,              ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_task_class,              ONLY: t_jsb_process_task, t_jsb_task_options
  USE mo_jsb_process_class,           ONLY: t_jsb_process, A2L_, SPQ_, Q_RAD_

  ! "Integrate" routines
  USE mo_spq_util,                    ONLY: reset_spq_fluxes_real => reset_spq_fluxes
  USE mo_q_spq_update_turbulence,     ONLY: update_soil_turbulence_real => update_soil_turbulence
  USE mo_q_spq_update_physics,        ONLY: update_spq_physics_real => update_spq_physics

  dsl4jsb_Use_memory(SPQ_)
  dsl4jsb_Use_memory(Q_RAD_)
  dsl4jsb_Use_memory(A2L_)

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Register_spq_tasks_quincy

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: reset_spq_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_reset_spq_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => reset_spq_fluxes      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_reset_spq_fluxes !< Aggregates computed task variables
  END TYPE tsk_reset_spq_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_soil_turbulence task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_turbulence
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_turbulence      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_turbulence   !< Aggregates computed task variables
  END TYPE tsk_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  !> Type definition: update_spq_physics task
  !!
  !-----------------------------------------------------------------------------------------------------
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_spq_physics
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_spq_physics      !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_spq_physics   !< Aggregates computed task variables
  END TYPE tsk_spq_physics


  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: reset_spq_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_reset_spq_fluxes
    PROCEDURE Create_task_reset_spq_fluxes         !< Constructor function for task
  END INTERFACE tsk_reset_spq_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_soil_turbulence task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_soil_turbulence
    PROCEDURE Create_task_update_soil_turbulence         !< Constructor function for task
  END INTERFACE tsk_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  !> Constructor interface: update_spq_physics task
  !!
  !-----------------------------------------------------------------------------------------------------
  INTERFACE tsk_spq_physics
    PROCEDURE Create_task_update_spq_physics         !< Constructor function for task
  END INTERFACE tsk_spq_physics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_spq_interface'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Register tasks: SPQ_ process
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE Register_spq_tasks_quincy(this, model_id)
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_reset_spq_fluxes(model_id))
    CALL this%Register_task(tsk_soil_turbulence(model_id))
    CALL this%Register_task(tsk_spq_physics(model_id))
  END SUBROUTINE Register_spq_tasks_quincy

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: reset_spq_fluxes task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_reset_spq_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_reset_spq_fluxes::return_ptr)
    CALL return_ptr%Construct(name='reset_spq_fluxes', process_id=SPQ_, owner_model_id=model_id)

  END FUNCTION Create_task_reset_spq_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_soil_turbulence task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_soil_turbulence(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_turbulence::return_ptr)
    CALL return_ptr%Construct(name='soil_turbulence', process_id=SPQ_, owner_model_id=model_id)

  END FUNCTION Create_task_update_soil_turbulence

  !-----------------------------------------------------------------------------------------------------
  !> Constructor: update_spq_physics task
  !!
  !-----------------------------------------------------------------------------------------------------
  FUNCTION Create_task_update_spq_physics(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_spq_physics::return_ptr)
    CALL return_ptr%Construct(name='spq_physics', process_id=SPQ_, owner_model_id=model_id)

  END FUNCTION Create_task_update_spq_physics

  ! ------------------------------------------------------------------------------------------------------- !
  ! Wrapper for update routine from different module
  ! ------------------------------------------------------------------------------------------------------- !
  SUBROUTINE update_spq_physics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_spq_physics_real(tile, options)

  END SUBROUTINE update_spq_physics

  SUBROUTINE update_soil_turbulence(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL update_soil_turbulence_real(tile, options)

  END SUBROUTINE update_soil_turbulence

  SUBROUTINE reset_spq_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CALL reset_spq_fluxes_real(tile, options)

  END SUBROUTINE reset_spq_fluxes

  ! ======================================================================================================= !
  !>Implementation of "aggregate": reset_spq_fluxes task
  !>
  SUBROUTINE aggregate_reset_spq_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SPQ_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_reset_spq_fluxes'
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
    dsl4jsb_Get_memory(SPQ_)

    ! nothing to aggregate after this routine

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_reset_spq_fluxes

  ! ======================================================================================================= !
  !>Implementation of "aggregate": soil_turbulence task
  !>
  SUBROUTINE aggregate_soil_turbulence(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_turbulence'
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
    dsl4jsb_Get_memory(A2L_)

    ! dsl4jsb_Aggregate_onChunk(A2L_, drag_srf             , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(A2L_, pch                  , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(A2L_, t_acoef              , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(A2L_, t_bcoef              , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(A2L_, q_acoef              , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(A2L_, q_bcoef              , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_turbulence

  ! ======================================================================================================= !
  !>Implementation of "aggregate": spq_physics task
  !>
  SUBROUTINE aggregate_spq_physics(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER                :: model
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(SPQ_)
    ! dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_RAD_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_aggregator),  POINTER         :: weighted_by_fract
    INTEGER                                   :: iblk, ics, ice
    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_spq_physics'
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
    dsl4jsb_Get_memory(SPQ_)
    ! dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_RAD_)

    ! dsl4jsb_Aggregate_onChunk(A2L_, pch                   , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(A2L_, drag_srf              , weighted_by_fract)
    ! ---------------------------
    ! dsl4jsb_Aggregate_onChunk(L2A_, fact_q_air            , weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(L2A_, fact_qsat_srf         , weighted_by_fract)
    ! ---------------------------
    dsl4jsb_Aggregate_onChunk(Q_RAD_, net_radiation           , weighted_by_fract)
    ! ---------------------------
    ! 2D
    dsl4jsb_Aggregate_onChunk(SPQ_, fact_q_air              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, fact_qsat_srf           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, t_srf_old               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, zril_old                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, transpiration           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_soil_root_theta       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_soil_root_pot         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_skin                  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_soil_root             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, t_soil_root             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, pet                     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, interception            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, evaporation             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, evaporation_snow        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, srf_runoff              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, drainage                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, drainage_fraction       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, ground_heat_flx_old     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, ground_heat_flx         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, latent_heat_flx         , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, sensible_heat_flx       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, snow_height             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, snow_soil_heat_flux     , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, snow_melt_to_soil       , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, snow_srf_heat_flux      , weighted_by_fract)
    ! 3D
    dsl4jsb_Aggregate_onChunk(SPQ_, w_soil_sl               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_soil_pot_sl           , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, gw_runoff_sl            , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, percolation_sl          , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, frac_w_lat_loss_sl      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, drainage_sl             , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, t_soil_sl               , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_ice_sl                , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_soil_freeze_flux      , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_soil_melt_flux        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_snow_snl              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, t_snow_snl              , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, snow_lay_thickness_snl  , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, snow_present_snl        , weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SPQ_, w_liquid_snl            , weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_spq_physics

#endif
END MODULE mo_q_spq_interface
