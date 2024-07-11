!> Contains the interfaces to the assimi process
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
!NEC$ options "-finline-file=externals/jsbach/src/assimilation/mo_assimi_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/hydrology/mo_hydro_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_assimi_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,     ONLY: debug_on
  USE mo_jsb_time,        ONLY: is_time_experiment_start
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_grid_class,     ONLY: t_jsb_grid
  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  !USE mo_jsb_config_class,   ONLY: t_jsb_config, t_jsb_config_p
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module (Get_assimi_memory and Get_assimi_config)
  dsl4jsb_Use_processes ASSIMI_, HYDRO_, PHENO_, RAD_, A2L_, SEB_

  ! Use of process configurations (t_assimi_config)
  dsl4jsb_Use_config(ASSIMI_)

  ! Use of process memories (t_assimi_memory)
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(ASSIMI_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(PHENO_)
  dsl4jsb_Use_memory(RAD_)
  dsl4jsb_Use_memory(SEB_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  Register_assimi_tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_interface'

  !> Type definition for assimilation_scaling_factors task
  TYPE, EXTENDS(t_jsb_process_task) ::   tsk_assimilation_scaling_factors
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_assimilation_scaling_factors    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimilation_scaling_factors !< Aggregates computed task variables
  END TYPE   tsk_assimilation_scaling_factors

  !> Constructor interface for assimilation_scaling_factors task
  INTERFACE   tsk_assimilation_scaling_factors
    PROCEDURE Create_task_assimilation_scaling_factors         !< Constructor function for task
  END INTERFACE   tsk_assimilation_scaling_factors

  !> Type definition for canopy_cond_unstressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_unstressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_unstressed    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_unstressed !< Aggregates computed task variables
  END TYPE tsk_canopy_cond_unstressed

  !> Constructor interface for canopy_cond_unstressed task
  INTERFACE tsk_canopy_cond_unstressed
    PROCEDURE Create_task_canopy_cond_unstressed     !< Constructor function for task
  END INTERFACE tsk_canopy_cond_unstressed

  !> Type definition for canopy_cond_stressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_stressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_stressed    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_stressed !< Aggregates computed task variables
  END TYPE tsk_canopy_cond_stressed

  !> Constructor interface for canopy_cond_stressed task
  INTERFACE tsk_canopy_cond_stressed
    PROCEDURE Create_task_canopy_cond_stressed     !< Constructor function for task
  END INTERFACE tsk_canopy_cond_stressed

  !> Type definition for assimilation task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_assimilation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_assimilation    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_assimilation !< Aggregates computed task variables
  END TYPE tsk_assimilation

  !> Constructor interface for assimilation task
  INTERFACE tsk_assimilation
    PROCEDURE Create_task_assimilation     !< Constructor function for task
  END INTERFACE tsk_assimilation

  !> Type definition for NPP buffer task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_npp_buffer
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_npp_buffer    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_npp_buffer !< Aggregates computed task variables
  END TYPE tsk_npp_buffer

  !> Constructor interface for assimilation task
  INTERFACE tsk_npp_buffer
    PROCEDURE Create_task_npp_buffer     !< Constructor function for task
  END INTERFACE tsk_npp_buffer

CONTAINS

  ! ------------------------------------------------------------------------------------------------------- !
  !> Constructors for tasks
  !>

  ! ======================================================================================================= !
  !> Constructor for assimilation_scaling_factors task
  !
  ! @param[in]     model_id     Model id
  ! @return        return_ptr   Instance of process task "assimi"
  FUNCTION Create_task_assimilation_scaling_factors(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimilation_scaling_factors::return_ptr)
    CALL return_ptr%Construct(name='assimilation_scaling_factors', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_assimilation_scaling_factors

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_unstressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_unstressed"
  !!
  FUNCTION Create_task_canopy_cond_unstressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_unstressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_unstressed', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for canopy_cond_stressed task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "canopy_cond_stressed"
  !!
  FUNCTION Create_task_canopy_cond_stressed(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_canopy_cond_stressed::return_ptr)
    CALL return_ptr%Construct(name='canopy_cond_stressed', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for assimilation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "assimilation"
  !!
  FUNCTION Create_task_assimilation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_assimilation::return_ptr)
    CALL return_ptr%Construct(name='assimilation', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_assimilation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for npp_buffer task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "npp_buffer"
  !!
  FUNCTION Create_task_npp_buffer(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_npp_buffer::return_ptr)
    CALL return_ptr%Construct(name='npp_buffer', process_id=ASSIMI_, owner_model_id=model_id)

  END FUNCTION Create_task_npp_buffer

  ! ======================================================================================================= !
  !> Register tasks for assimi process
  !>
  ! TODO docu:
  !
  ! @param[in,out] this      Instance of assimi process class
  ! @param[in]     model_id  Model id
  SUBROUTINE Register_assimi_tasks(this, model_id)
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model

    ! in/out
    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,                 INTENT(in) :: model_id

    CALL this%Register_task(tsk_assimilation_scaling_factors(model_id))
    CALL this%Register_task(tsk_canopy_cond_unstressed(model_id))
    CALL this%Register_task(tsk_canopy_cond_stressed(model_id))
    CALL this%Register_task(tsk_assimilation(model_id))
    CALL this%Register_task(tsk_npp_buffer(model_id))
  END SUBROUTINE Register_assimi_tasks

  ! ================================================================================================================================
  !>
  !> Implementation to calculate the height scaling factors of plant enzyms for assimilation
  !! Task "assimi" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  ! The scaling factors represent the enzymatic distribution of plants in their canopy layers. This correlates
  ! with the amount of nitrogen in these layers (therefore called nitrogen_scaling_factors in JSBACH3).
  ! Eqs. (106), (107) in Knorr (1997)
  ! Note that this routine should only be applied to special vegetation classes
  SUBROUTINE update_assimilation_scaling_factors(tile, options)

    USE mo_orbit_solar,    ONLY: inquire_declination
    USE mo_jsb_grid,       ONLY: Get_grid
    USE mo_assimi_process, ONLY: calc_assimilation_scaling_factors

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER   :: model
    TYPE(t_jsb_grid),  POINTER   :: grid
    INTEGER                      :: iblk, ics, ice, nc, ic
    REAL(wp),          POINTER   :: lat(:)
    REAL(wp)                     :: declination ! Solar Declination angle
    INTEGER                      :: i
    INTEGER                      :: ncanopy          ! number of canopy layers
    REAL(wp),          POINTER   :: canopy_bound_lai(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':assimilation_scaling_factors'

    ! Declare process configuration and memory Pointers
    dsl4jsb_Def_config(ASSIMI_)
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers to variables in memory
    dsl4jsb_Real3D_onChunk :: scaling_fact_cl
    dsl4jsb_Real2D_onChunk :: lai

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get process config
    dsl4jsb_Get_config(ASSIMI_)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(PHENO_)

    ! Get process variables
    ncanopy          =  dsl4jsb_Config(ASSIMI_)%ncanopy
    canopy_bound_lai(0:ncanopy) => dsl4jsb_Config(ASSIMI_)%canopy_bound_lai(0:ncanopy)

    grid => Get_grid(model%grid_id)
    lat  => grid%lat(ics:ice, iblk)

    dsl4jsb_Get_var2D_onChunk(PHENO_,  lai)             ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_, scaling_fact_cl) ! out

    ! ---------------------------
    ! Go

    CALL inquire_declination(declination) ! out

    IF (dsl4jsb_Lctlib_param(NitrogenScalingFlag)) THEN
      DO i=1,ncanopy
        CALL calc_assimilation_scaling_factors( &
          & declination,                        & ! input
          & canopy_bound_lai(i-1),              & ! input
          & lai(:),                             & ! input
          & lat(:),                             & ! input
          & scaling_fact_cl(:,i)                & ! output
          & )
      END DO
    ELSE                                           ! to some vegetation classes nitrogen scaling is not applied
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO i=1,ncanopy
        !$ACC LOOP GANG VECTOR
        DO ic=1,nc
          scaling_fact_cl(ic,i) = 1._wp
        END DO
      END DO
      !$ACC END PARALLEL
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_assimilation_scaling_factors

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "assimilation_scaling_factors"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_assimilation_scaling_factors(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_assimilation_scaling_factors'

    INTEGER  :: iblk, ics, ice !, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, scaling_fact_cl, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_assimilation_scaling_factors

  ! ================================================================================================================================
  !>
  !> Implementation of task canopy_cond_unstressed
  !! Task "canopy_cond_unstressed" calculates the carbon assimilation of plants and canopy conductance without water stress.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_unstressed(tile, options)

    USE mo_assimi_process,      ONLY: calc_assimilation_waterunlimited

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER                :: model
    !REAL(wp)                                  :: dtime
    INTEGER                                   :: iblk, ics, ice, nc, ic
    INTEGER                                   :: icanopy

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_unstressed'

    ! Declare process memories
    dsl4jsb_Def_memory(ASSIMI_)
    !dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(A2L_)

    dsl4jsb_Def_config(ASSIMI_)

    ! Declare pointers to variables in memory
    INTEGER :: ncanopy              ! number of canopy layers
    REAL(wp) :: &
      & CarboxRate, &
      & ETransport
    LOGICAL :: C4Flag

    dsl4jsb_Real2D_onChunk :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: press_srf             ! Surface pressure
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.
    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited ! Canopy conductance no limitation
    !dsl4jsb_Real2D_onChunk :: canopy_cond_limited   ! Canopy conductance water limited

    dsl4jsb_Real2D_onChunk :: CO2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    dsl4jsb_Real3D_onChunk :: scaling_fact_cl
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl
    dsl4jsb_Real3D_onChunk :: lai_cl                  ! Leaf area index [-]


    ! Get local variables from options argument
    !dtime   = options%dtime
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(ASSIMI_)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    !dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(A2L_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_,      t_air)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)               ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      CO2_air_mol)             ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      par_down_mol)            ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      apar_per_lai_cl)         ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      lai_cl)                  ! in

    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   scaling_fact_cl)         ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   canopy_cond_cl)          ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   canopy_cond_unlimited)   ! out

    ncanopy     = dsl4jsb_Config(ASSIMI_)%ncanopy
    C4Flag      = dsl4jsb_Lctlib_param(C4Flag)
    CarboxRate  = dsl4jsb_Lctlib_param(CarboxRate)
    ETransport  = dsl4jsb_Lctlib_param(ETransport)

    ! TODO: LOOP SEQ over icanopy would create "complex loop carried dependence ...", why?
    DO icanopy=1,ncanopy
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        canopy_cond_cl(ic,icanopy) = calc_assimilation_waterunlimited( &
        & C4Flag,                      &
        & CarboxRate,                  &
        & ETransport,                  &
        & t_air(ic),                   &
        & press_srf(ic),               &
        & par_down_mol(ic),            &
        & apar_per_lai_cl(ic,icanopy), &
        & CO2_air_mol(ic),             &
        & scaling_fact_cl(ic,icanopy)  &
        & )
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL
    END DO

    ! Compute unlimited conductance for the whole canopy (all layers) considering the lai
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      canopy_cond_unlimited(ic) = MAX(1.e-20_wp,SUM(canopy_cond_cl(ic,:) * lai_cl(ic,:)))
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_unstressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_unstressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)
    !dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_unstressed'

    INTEGER  :: iblk, ics, ice !, nc

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)
    !dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_unlimited, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "update" for task "canopy_cond_stressed"
  !! Task "update_canopy_cond_stressed" calculates canopy conductance under water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_stressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance
    USE mo_phy_schemes,   ONLY: qsat_water

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(ASSIMI_)

    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf

    INTEGER :: iblk, ics, ice, nc, ic
    INTEGER :: icanopy, ncanopy
    REAL(wp) :: ratio_lim_over_unlim(options%nc)
    LOGICAL :: q_air_gt_qsat_tmp
    LOGICAL :: use_tmx

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(ASSIMI_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    use_tmx = model%config%use_tmx

    dsl4jsb_Get_config(ASSIMI_)

    ncanopy = dsl4jsb_Config(ASSIMI_)%ncanopy

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      q_air)               ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)           ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)        ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   canopy_cond_unlimited) ! in
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   canopy_cond_limited) ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   canopy_cond_cl)       ! inout

    ! Compute (actual) canopy (stomatal) conductance under water stress.
    !$ACC DATA CREATE(ratio_lim_over_unlim)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
    !$ACC   PRIVATE(q_air_gt_qsat_tmp)
    DO ic = 1, nc
      q_air_gt_qsat_tmp = q_air(ic) > qsat_water(t(ic),press_srf(ic), use_convect_tables=.NOT. use_tmx)
      canopy_cond_limited(ic) = &
        &  get_canopy_conductance(canopy_cond_unlimited(ic), & ! in, unstressed canopy conductance
                                  water_stress(ic),          & ! in, water stress factor
                                  q_air_gt_qsat_tmp          & ! in, atmosphere saturated?
                                 )
      ratio_lim_over_unlim(ic) = canopy_cond_limited(ic) / canopy_cond_unlimited(ic)
    END DO
    !$ACC END PARALLEL LOOP

    ! Scale unlimited conductance per leaf area with the ratio of limited to unlimited conductance for
    ! whole canopy, ie. compute new limited conductance per leaf area and canopy layer
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    !$ACC LOOP SEQ
    DO icanopy=1,ncanopy
      !$ACC LOOP VECTOR
      DO ic=1,nc
        ! canopy_cond_cl(ic,icanopy) = canopy_cond_cl(ic,icanopy) * canopy_cond_limited(ic) / canopy_cond_unlimited(ic)
         canopy_cond_cl(ic,icanopy) = canopy_cond_cl(ic,icanopy) * ratio_lim_over_unlim(ic)
      END DO
      !$ACC END LOOP
    END DO
    !$ACC END LOOP
    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_stressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_stressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_stressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_limited, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl,      weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_stressed

  ! ================================================================================================================================
  !>
  !> Implementation of assimilation for water-limited conditions
  !! Task "assimilation" calculates the carbon assimilation of plants.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_assimilation(tile, options)

    USE mo_assimi_process,  ONLY: calc_assimilation_waterlimited, calc_NPP_pot_rate

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER                    :: iblk, ics, ice, nc, ic
    REAL(wp)                   :: dtime
    INTEGER                    :: icanopy

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_assimilation'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_memory(RAD_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_config(ASSIMI_)

    INTEGER :: ncanopy                 ! number of canopy layers

    ! Declare pointers to variables in memory

    dsl4jsb_Real2D_onChunk :: t_air                 ! Atmosphere temperature (lowest layer) in Kelvin!
    dsl4jsb_Real2D_onChunk :: press_srf             ! Surface pressure
    dsl4jsb_Real2D_onChunk :: par_down_mol          ! Downward PAR flux in mol (photons)/(m^2 s)
    dsl4jsb_Real3D_onChunk :: apar_per_lai_cl       ! Input. Absorbed PAR of canopy layer.

    dsl4jsb_Real2D_onChunk :: CO2_air_mol           ! CO2 particle mixing ratio [molCO2/molDryAir]

    dsl4jsb_Real2D_onChunk :: veg_fract_correction  ! factor to convert from tile to canopy area
    dsl4jsb_Real2D_onChunk :: gross_assimilation_ca
    dsl4jsb_Real2D_onChunk :: gross_assimilation
    dsl4jsb_Real2D_onChunk :: dark_respiration_ca
    dsl4jsb_Real3D_onChunk :: dark_respiration_cl
    dsl4jsb_Real2D_onChunk :: dark_respiration
    dsl4jsb_Real2D_onChunk :: NPP_pot_rate_ca
    dsl4jsb_Real3D_onChunk :: scaling_fact_cl        ! Input
    dsl4jsb_Real2D_onChunk :: fract_fpc_max          ! in
    dsl4jsb_Real3D_onChunk :: canopy_cond_cl         ! InOut for water limited case here.
    dsl4jsb_Real3D_onChunk :: gross_assimilation_cl  ! Output.
    dsl4jsb_Real3D_onChunk :: carbox_rate_max_cl     ! Output. Maximum carboxylation rate of canopy layer (=VCmax)
                                                         ! [Mol(CO2)/m^2(canopy)/s]. carbox_rate_max_leaf adapted to canopy layer
                                                         ! and actual temperature.
    dsl4jsb_Real3D_onChunk :: e_transport_rate_max_cl ! Output
    dsl4jsb_Real3D_onChunk :: carbox_rate_cl          ! Output
    dsl4jsb_Real3D_onChunk :: e_transport_rate_cl     ! Output
    dsl4jsb_Real3D_onChunk :: CO2_conc_leaf_cl        ! CO2 concentration inside leaf  [mol(CO2)/mol(Air)]:
                                                      ! Output only for the water limited case here
    dsl4jsb_Real3D_onChunk :: lai_cl                  ! Leaf area index [-]
    dsl4jsb_Real2D_onChunk :: day_NPP_sum

    ! Set local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set process config and memories
    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ASSIMI_)

    ncanopy = dsl4jsb_Config(ASSIMI_)%ncanopy

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)
    dsl4jsb_Get_memory(RAD_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(PHENO_)

    ! Get process variables
    dsl4jsb_Get_var2D_onChunk(A2L_,      t_air)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      CO2_air_mol)             ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)               ! in
    dsl4jsb_Get_var2D_onChunk(RAD_,      par_down_mol)            ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      apar_per_lai_cl)         ! in
    dsl4jsb_Get_var3D_onChunk(RAD_,      lai_cl)                  ! in

    dsl4jsb_Get_var2D_onChunk(PHENO_,    veg_fract_correction)    ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   scaling_fact_cl)         ! in
    dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc_max )          ! in
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   canopy_cond_cl)          ! in

    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   gross_assimilation_ca)   ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   gross_assimilation)      ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   dark_respiration_ca)     ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   dark_respiration)        ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   NPP_pot_rate_ca)         ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   gross_assimilation_cl)   ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   dark_respiration_cl)     ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   carbox_rate_max_cl)      ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   e_transport_rate_max_cl) ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   e_transport_rate_cl)     ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   carbox_rate_cl)          ! out
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,   CO2_conc_leaf_cl)        ! out

    dsl4jsb_Get_var2D_onChunk(PHENO_,    day_NPP_sum)             ! inout

    ! ---------------------------
    ! Go

    DO icanopy=1,ncanopy
      CALL calc_assimilation_waterlimited(    &
        & dsl4jsb_Lctlib_param(C4Flag),       &
        & dsl4jsb_Lctlib_param(CarboxRate),   &
        & dsl4jsb_Lctlib_param(ETransport),   &
        & t_air(:),                           & ! in
        & press_srf(:),                       & ! in
        & par_down_mol(:),                    & ! in
        & apar_per_lai_cl(:,icanopy),         & ! in
        & CO2_air_mol(:),                     & ! in
        & scaling_fact_cl(:,icanopy),         & ! in
        & CO2_conc_leaf_cl(:,icanopy),        & ! out only for the water limited case here
        & canopy_cond_cl(:,icanopy),          & ! in  only for the water limited case here
        & gross_assimilation_cl(:,icanopy),   & ! out
        & dark_respiration_cl(:,icanopy),     & ! out
        & carbox_rate_max_cl(:,icanopy),      & ! out
        & e_transport_rate_max_cl(:,icanopy), & ! out
        & carbox_rate_cl(:,icanopy),          & ! out
        & e_transport_rate_cl(:,icanopy)      & ! out
        & )
    END DO

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      ! convert from leaf area of layer to canopy ground area
      gross_assimilation_ca(ic) = SUM(gross_assimilation_cl(ic,:) * lai_cl(ic,:))
      ! convert from canopy area to tile box area
      gross_assimilation(ic) = gross_assimilation_ca(ic) * veg_fract_correction(ic) * fract_fpc_max(ic)

      dark_respiration_ca(ic) = SUM(dark_respiration_cl(ic,:) * lai_cl(ic,:))
      ! convert from canopy area to tile box area
      dark_respiration(ic) =  dark_respiration_ca(ic)  * veg_fract_correction(ic) * fract_fpc_max(ic)
    END DO
    !$ACC END PARALLEL LOOP

    CALL calc_NPP_pot_rate(gross_assimilation_ca(:), dark_respiration_ca(:), NPP_pot_rate_ca(:))

    ! Accumulate NPP_pot_rate_ca. day_NPP_sum is divided by the day length and reset to zero in
    ! update_phenology_logrop at the beginning of each day.
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      day_NPP_sum(ic) = day_NPP_sum(ic) + NPP_pot_rate_ca(ic) * dtime
    END DO
    !$ACC END PARALLEL LOOP

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_assimilation

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "assimilation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_assimilation(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice !, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_assimilation'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, canopy_cond_cl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation_ca,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, gross_assimilation,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, dark_respiration_ca,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, dark_respiration,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, NPP_pot_rate_ca,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, CO2_conc_leaf_cl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, carbox_rate_cl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, carbox_rate_max_cl,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, e_transport_rate_cl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, e_transport_rate_max_cl,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_assimilation

  ! ================================================================================================================================
  !>
  !> Implementation of NPP buffer
  !! Task "npp_buffer" calculates multi-annual average of potential  NPP for use in the NLCC process
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_npp_buffer(tile, options)

    USE mo_assimi_process,      ONLY: calc_NPPbuf
    USE mo_jsb_time,            ONLY: is_newyear

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    INTEGER                    :: iblk, ics, ice !, nc
    REAL(wp)                   :: dtime
    LOGICAL                    :: lstart, init_running_means, init_npp, new_year

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_npp_buffer'

    ! Declare pointers to process configuration and memory
    dsl4jsb_Def_memory(ASSIMI_)
    dsl4jsb_Def_config(ASSIMI_)

    ! Declare pointers to variables in memory

    dsl4jsb_Real2D_onChunk :: NPP_pot_rate_ca
    dsl4jsb_Real2D_onChunk :: seconds_year
    dsl4jsb_Real2D_onChunk :: NPP_sum_year
    dsl4jsb_Real2D_onChunk :: NPP_mean_5year
    dsl4jsb_Real2D_onChunk :: land_cover_class
    dsl4jsb_Real2D_onChunk :: bclimit_min_cold_mmtemp
    dsl4jsb_Real2D_onChunk :: bclimit_max_cold_mmtemp
    dsl4jsb_Real2D_onChunk :: bclimit_max_warm_mmtemp
    dsl4jsb_Real2D_onChunk :: bclimit_min_temprange
    dsl4jsb_Real2D_onChunk :: bclimit_min_gdd
    dsl4jsb_Real2D_onChunk :: tau_c_woods

    ! Set local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    !nc      = options%nc
    dtime   = options%dtime

    ! If process is not to be calculated on this tile, do nothing
    IF (.NOT. tile%Is_process_calculated(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Set process config and memories
    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(ASSIMI_)

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_)

    ! Get process variables

    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   NPP_pot_rate_ca)         ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   seconds_year)            ! inout
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   NPP_sum_year)            ! inout
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   NPP_mean_5year)          ! inout
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   land_cover_class)        ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   bclimit_min_cold_mmtemp) ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   bclimit_max_cold_mmtemp) ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   bclimit_max_warm_mmtemp) ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   bclimit_min_temprange)   ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   bclimit_min_gdd)         ! out
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,   tau_c_woods)             ! out

    ! Accumulate NPP_sum_year and update NPP_mean_5year for NLCC
    lstart = is_time_experiment_start(options%current_datetime)
    new_year = is_newyear(options%current_datetime, dtime)

    init_running_means = dsl4jsb_Config(ASSIMI_)%init_running_means
    init_NPP = .TRUE.
    IF (MINVAL(NPP_mean_5year(:)) < -999._wp) init_NPP = .FALSE.
    CALL calc_NPPbuf(lstart, init_running_means, new_year, dtime, NPP_pot_rate_ca, &
                     seconds_year, NPP_sum_year, NPP_mean_5year)

    IF (.NOT. init_NPP .AND. debug_on() .AND. iblk==1) THEN
       IF (MINVAL(NPP_mean_5year(:)) >= -999._wp) CALL message(TRIM(routine), 'initialisation of NPP_mean_5year')
    ENDIF

    ! write lct information on variables to make them available on pft-tiles via function collect_var for NLCC
    land_cover_class(:) = REAL(dsl4jsb_Lctlib_param(LandcoverClass), wp)
    bclimit_min_cold_mmtemp(:) = dsl4jsb_Lctlib_param(bclimit_min_cold_mmtemp)
    bclimit_max_cold_mmtemp(:) = dsl4jsb_Lctlib_param(bclimit_max_cold_mmtemp)
    bclimit_max_warm_mmtemp(:) = dsl4jsb_Lctlib_param(bclimit_max_warm_mmtemp)
    bclimit_min_temprange(:) = dsl4jsb_Lctlib_param(bclimit_min_temprange)
    bclimit_min_gdd(:) = dsl4jsb_Lctlib_param(bclimit_min_gdd)
    tau_c_woods(:) = dsl4jsb_Lctlib_param(tau_c_woods)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_npp_buffer

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "npp_buffer"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_npp_buffer(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(ASSIMI_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    ! Local variables
    INTEGER  :: iblk, ics, ice !, nc

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_npp_buffer'

    ! Get local variables from options argument
    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    !nc   = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(ASSIMI_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(ASSIMI_, seconds_year,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, NPP_sum_year,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, NPP_mean_5year,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, land_cover_class,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, bclimit_min_cold_mmtemp,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, bclimit_max_cold_mmtemp,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, bclimit_max_warm_mmtemp,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, bclimit_min_temprange,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, bclimit_min_gdd,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(ASSIMI_, tau_c_woods,              weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_npp_buffer

#endif
END MODULE mo_assimi_interface
