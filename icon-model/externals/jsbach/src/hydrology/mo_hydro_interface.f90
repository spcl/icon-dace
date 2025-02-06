!> Contains the interfaces to the hydro processes
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
!NEC$ options "-finline-file=externals/jsbach/src/hydrology/mo_hydro_process.pp-jsb.f90"
!NEC$ options "-finline-file=externals/jsbach/src/shared/mo_phy_schemes.pp-jsb.f90"
!NEC$ options "-finline-max-function-size=100"

MODULE mo_hydro_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  ! Use of basic structures
  USE mo_jsb_control,        ONLY: debug_on
  USE mo_jsb_time,           ONLY: is_time_experiment_start
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: message, message_text, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_grid, t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_grid, Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  ! Use of processes in this module
  dsl4jsb_Use_processes SEB_, SSE_, TURB_, PHENO_, A2L_, HYDRO_

  ! Use of process configurations
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)

  ! Use of process memories
  dsl4jsb_Use_memory(A2L_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(TURB_)
  dsl4jsb_Use_memory(PHENO_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_hydro_tasks, global_hydrology_diagnostics

  CHARACTER(len=*), PARAMETER :: modname = 'mo_hydro_interface'

  !> Type definition for surface_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_hydrology
  END TYPE tsk_surface_hydrology

  !> Constructor interface for surface_hydrology
  INTERFACE tsk_surface_hydrology
    PROCEDURE Create_task_surface_hydrology
  END INTERFACE tsk_surface_hydrology

  !> Type definition for soil_properties
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_properties
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_properties
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_properties
  END TYPE tsk_soil_properties

  !> Constructor interface for soil_properties
  INTERFACE tsk_soil_properties
    PROCEDURE Create_task_soil_properties
  END INTERFACE tsk_soil_properties

  !> Type definition for soil_hydrology
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_hydrology
  END TYPE tsk_soil_hydrology

  !> Constructor interface for soil_hydrology
  INTERFACE tsk_soil_hydrology
    PROCEDURE Create_task_soil_hydrology
  END INTERFACE tsk_soil_hydrology

  !> Type definition for evaporation task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_evaporation
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_evaporation
    PROCEDURE, NOPASS :: Aggregate => aggregate_evaporation
  END TYPE tsk_evaporation

  !> Constructor interface for evaporation task
  INTERFACE tsk_evaporation
    PROCEDURE Create_task_evaporation
  END INTERFACE tsk_evaporation

  !> Type definition for canopy_cond_unstressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_unstressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_unstressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_unstressed
  END TYPE tsk_canopy_cond_unstressed

  !> Constructor interface for canopy_cond_unstressed task
  INTERFACE tsk_canopy_cond_unstressed
    PROCEDURE Create_task_canopy_cond_unstressed
  END INTERFACE tsk_canopy_cond_unstressed

  !> Type definition for water_stress task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_stress
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_stress
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_stress
  END TYPE tsk_water_stress

  !> Constructor interface for water_stress task
  INTERFACE tsk_water_stress
    PROCEDURE Create_task_water_stress
  END INTERFACE tsk_water_stress

  !> Type definition for canopy_cond_stressed task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_canopy_cond_stressed
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_canopy_cond_stressed
    PROCEDURE, NOPASS :: Aggregate => aggregate_canopy_cond_stressed
  END TYPE tsk_canopy_cond_stressed

  !> Constructor interface for canopy_cond_stressed task
  INTERFACE tsk_canopy_cond_stressed
    PROCEDURE Create_task_canopy_cond_stressed
  END INTERFACE tsk_canopy_cond_stressed

  !> Type definition for snow_hydrology task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_ice_hydrology
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_ice_hydrology
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_ice_hydrology
  END TYPE tsk_snow_and_ice_hydrology

  !> Constructor interface for snow_hydrology task
  INTERFACE tsk_snow_and_ice_hydrology
    PROCEDURE Create_task_snow_and_ice_hydrology
  END INTERFACE tsk_snow_and_ice_hydrology

  !> Type definition for snow_and_wet_fraction task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snow_and_wet_fraction
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snow_and_wet_fraction    !< Advances task computation for one timestep
    PROCEDURE, NOPASS :: Aggregate => aggregate_snow_and_wet_fraction !< Aggregates computed task variables
  END TYPE tsk_snow_and_wet_fraction

  !> Constructor interface for snow_and_wet_fraction task
  INTERFACE tsk_snow_and_wet_fraction
    PROCEDURE Create_task_snow_and_wet_fraction                       !< Constructor function for task
  END INTERFACE tsk_snow_and_wet_fraction

  !> Type definition for water_balance task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_water_balance
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_water_balance
    PROCEDURE, NOPASS :: Aggregate => aggregate_water_balance
  END TYPE tsk_water_balance

  !> Constructor interface for water balance task
  INTERFACE tsk_water_balance
    PROCEDURE Create_task_water_balance
  END INTERFACE tsk_water_balance

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_hydrology"
  !!
  FUNCTION Create_task_surface_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_hydrology::return_ptr)
    CALL return_ptr%Construct(name='surface_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_properties task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_properties"
  !!
  FUNCTION Create_task_soil_properties(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_properties::return_ptr)
    CALL return_ptr%Construct(name='soil_properties', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_properties

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_hydrology"
  !!
  FUNCTION Create_task_soil_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_hydrology::return_ptr)
    CALL return_ptr%Construct(name='soil_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_hydrology

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
    CALL return_ptr%Construct(name='canopy_cond_unstressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_stress task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_stress"
  !!
  FUNCTION Create_task_water_stress(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_stress::return_ptr)
    CALL return_ptr%Construct(name='water_stress', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_stress

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
    CALL return_ptr%Construct(name='canopy_cond_stressed', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for evaporation task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "evaporation"
  !!
  FUNCTION Create_task_evaporation(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_evaporation::return_ptr)
    CALL return_ptr%Construct(name='evaporation', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_evaporation

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow hydrology task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "snow_hydrology"
  !!
  FUNCTION Create_task_snow_and_ice_hydrology(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snow_and_ice_hydrology::return_ptr)
    CALL return_ptr%Construct(name='snow_ice_hydrology', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_ice_hydrology

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snow_and_wet_fraction task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "tsk_snow_and_wet_fraction"
  !!
  FUNCTION Create_task_snow_and_wet_fraction(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snow_and_wet_fraction::return_ptr)
    CALL return_ptr%Construct(name='snow_and_wet_fraction', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_snow_and_wet_fraction

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for water_balance task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "water_balance"
  !!
  FUNCTION Create_task_water_balance(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_water_balance::return_ptr)
    CALL return_ptr%Construct(name='water_balance', process_id=HYDRO_, owner_model_id=model_id)

  END FUNCTION Create_task_water_balance

  ! -------------------------------------------------------------------------------------------------------
  !> Register tasks for hydrology process
  !!
  !! @param[in]     model_id  Model id
  !! @param[in,out] this      Instance of HYDRO_ process class
  !!
  SUBROUTINE Register_hydro_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task( tsk_surface_hydrology      (model_id))
    CALL process%Register_task( tsk_soil_properties        (model_id))
    CALL process%Register_task( tsk_soil_hydrology         (model_id))
    CALL process%Register_task( tsk_water_stress           (model_id))
    CALL process%Register_task( tsk_canopy_cond_unstressed (model_id))
    CALL process%Register_task( tsk_canopy_cond_stressed   (model_id))
    CALL process%Register_task( tsk_evaporation            (model_id))
    !CALL process%Register_task( tsk_snow_and_ice_hydrology (model_id))
    CALL process%Register_task( tsk_snow_and_wet_fraction  (model_id))
    CALL process%Register_task( tsk_water_balance          (model_id))

  END SUBROUTINE Register_hydro_tasks
  !
  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "surface_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_hydrology(tile, options)

    USE mo_hydro_process,          ONLY: get_soilhyd_properties, calc_surface_hydrology_land, &
      &                                  calc_surface_hydrology_glacier
    USE mo_jsb_physical_constants, ONLY: dens_snow
    USE mo_hydro_util,             ONLY: get_amount_in_rootzone

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(A2L_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: q_snocpymlt
    dsl4jsb_Real2D_onChunk :: t_air
    dsl4jsb_Real2D_onChunk :: wind_10m
    dsl4jsb_Real2D_onChunk :: rain
    dsl4jsb_Real2D_onChunk :: snow
    dsl4jsb_Real2D_onChunk :: wtr_skin
    dsl4jsb_Real2D_onChunk :: fract_skin
    dsl4jsb_Real2D_onChunk :: weq_snow_soil
    dsl4jsb_Real2D_onChunk :: snow_soil_dens
    dsl4jsb_Real2D_onChunk :: steepness
    dsl4jsb_Real2D_onChunk :: fract_pond
    dsl4jsb_Real2D_onChunk :: weq_pond_max
    dsl4jsb_Real2D_onChunk :: wtr_pond
    dsl4jsb_Real2D_onChunk :: ice_pond
    dsl4jsb_Real2D_onChunk :: weq_pond
    dsl4jsb_Real2D_onChunk :: pond_freeze
    dsl4jsb_Real2D_onChunk :: pond_melt
    dsl4jsb_Real2D_onChunk :: wtr_pond_net_flx
    dsl4jsb_Real2D_onChunk :: weq_snow
    dsl4jsb_Real2D_onChunk :: snow_accum
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: weq_snow_can
    dsl4jsb_Real2D_onChunk :: evapotrans_soil
    dsl4jsb_Real2D_onChunk :: evapo_skin
    dsl4jsb_Real2D_onChunk :: evapo_deficit
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: evapotrans
    dsl4jsb_Real2D_onChunk :: transpiration
    dsl4jsb_Real2D_onChunk :: evapopot
    dsl4jsb_Real2D_onChunk :: evapo_snow
    dsl4jsb_Real2D_onChunk :: evapo_pond
    dsl4jsb_Real2D_onChunk :: snowmelt
    dsl4jsb_Real2D_onChunk :: infilt
    dsl4jsb_Real2D_onChunk :: runoff
    dsl4jsb_Real2D_onChunk :: runoff_horton
    dsl4jsb_Real2D_onChunk :: drainage
    dsl4jsb_Real2D_onChunk :: runoff_glac
    dsl4jsb_Real2D_onChunk :: water_to_soil
    dsl4jsb_Real2D_onChunk :: weq_glac
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl
    dsl4jsb_Real3D_onChunk :: vol_wres_sl
    dsl4jsb_Real3D_onChunk :: hyd_cond_sat_sl
    dsl4jsb_Real3D_onChunk :: matric_pot_sl
    dsl4jsb_Real3D_onChunk :: bclapp_sl
    dsl4jsb_Real3D_onChunk :: pore_size_index_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: t_soil_sl

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & skinres_max,        &
      & skinres_canopy_max, &
      & weq_rootzone,       &
      & weq_rootzone_max,   &
      & weq_soil,           &
      & trans_tmp
    REAL(wp), DIMENSION(options%nc, options%nsoil_w) :: &
      & weq_soil_sl_tmp,           &
      & weq_wsat_soil_sl_tmp,      &
      & weq_wres_soil_sl_tmp,      &
      & ice_impedance_sl

    INTEGER  :: iblk, ics, ice, nc, ic, nsoil, is
    REAL(wp) :: dtime
    ! Variables from configs
    REAL(wp) :: &
      & w_skin_max_config, &
      & snow_depth_max_config
    INTEGER  :: &
      & hydro_scale
    LOGICAL  :: &
      & is_experiment_start,    &
      & l_compat401_config,     &
      & l_dynsnow_config,       &
      & l_infil_subzero_config, &
      & ltpe_closed,            &
      & use_tmx

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    nsoil = options%nsoil_w
    dtime = options%dtime
    is_experiment_start = is_time_experiment_start(options%current_datetime)

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    l_compat401_config = model%config%l_compat401

    use_tmx = model%config%use_tmx

    dsl4jsb_Get_config(HYDRO_)
    w_skin_max_config      = dsl4jsb_Config(HYDRO_)%w_skin_max ! Maximum capacity of skin reservoir (soil + canopy)
    snow_depth_max_config  = dsl4jsb_Config(HYDRO_)%snow_depth_max
    l_infil_subzero_config = dsl4jsb_Config(HYDRO_)%l_infil_subzero
    hydro_scale            = dsl4jsb_Config(HYDRO_)%hydro_scale ! choice of hydrological scale scheme

    dsl4jsb_Get_config(SSE_)
    l_dynsnow_config = dsl4jsb_Config(SSE_)%l_dynsnow

    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)

    ltpe_closed = .FALSE.
    IF (model%config%tpe_scheme == 'closed') ltpe_closed = .TRUE.

    IF (tile%is_lake) THEN
      ! compute runoff (P-E) for lakes here and exit
      ! Note: for the purpose of the water balance it is assumed here that all rain and snow fall is going into runoff
      ! @todo: compute actual water balance (snow/ice melt) in lake model and use for runoff. The corresponding energy
      ! fluxes and snow/ice budgets are already computed in the lake model.
      dsl4jsb_Get_var2D_onChunk(A2L_,   rain)               ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,   snow)               ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)           ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_, infilt)             ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)             ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)           ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)           ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_snow)         ! out
      IF (use_tmx) dsl4jsb_Get_var2D_onChunk(HYDRO_, q_snocpymlt)    ! out

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
      DO ic=1,nc
        runoff  (ic) = rain(ic) + snow(ic) + evapotrans(ic)
        infilt(ic)   = 0._wp
        drainage(ic) = 0._wp
        weq_snow(ic) = 0._wp
        evapo_snow(ic) = evapopot(ic)
      END DO
      !$ACC END PARALLEL LOOP
      IF (use_tmx) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) ASYNC(1) GANG VECTOR
        DO ic=1,nc
          q_snocpymlt(ic) = 0._wp
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      !$ACC WAIT(1)

      RETURN
    END IF

    dsl4jsb_Get_memory(SSE_)
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    dsl4jsb_Get_var2D_onChunk(A2L_,   t_air)          ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   wind_10m)       ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   rain)           ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,   snow)           ! in

    dsl4jsb_Get_var2D_onChunk(SSE_,   hcap_grnd)      ! in

    dsl4jsb_Get_var2D_onChunk(HYDRO_, steepness)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow)       ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens) ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapotrans)     ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapopot)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_snow)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, q_snocpymlt)    ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, infilt)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff)         ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, runoff_horton)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, drainage)       ! out
    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)       ! in

    IF (tile%is_glacier) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_glac)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   runoff_glac)      ! inout
    ELSE
      dsl4jsb_Get_var3D_onChunk(SSE_,     t_soil_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   wtr_soil_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   ice_soil_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   wtr_soil_pot_scool_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   vol_porosity_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   vol_wres_sl)      ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   hyd_cond_sat_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   matric_pot_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   bclapp_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   pore_size_index_sl)     ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,   soil_depth_sl)    ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   fract_skin)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_skin)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   water_to_soil)    ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   snow_accum)       ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapotrans_soil)  ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_skin)       ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_deficit)    ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   fract_pond)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_pond_max)     ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_pond)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   ice_pond)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_pond)         ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_pond)       ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   pond_freeze)      ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   pond_melt)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   wtr_pond_net_flx) ! out
      IF (tile%contains_vegetation) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,   transpiration)    ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,   weq_snow_can)     ! in
        dsl4jsb_Get_var2D_onChunk(PHENO_,   lai)              ! in
        dsl4jsb_Get_var2D_onChunk(PHENO_,   fract_fpc_max)    ! in
      END IF
    END IF

    IF (is_experiment_start) THEN ! Start of experiment
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        weq_snow(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (tile%is_glacier) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        CALL calc_surface_hydrology_glacier( &
          & is_experiment_start,             & ! in
          & dtime,                           & ! in
          & t_unfilt(ic),                    & ! in
          & fract_snow(ic),                  & ! in
          & hcap_grnd(ic),                   & ! in
          & evapotrans(ic),                  & ! in
          & evapopot(ic),                    & ! in
          & rain(ic),                        & ! in
          & snow(ic),                        & ! in
          & weq_glac(ic),                    & ! inout
          & q_snocpymlt(ic),                 & ! out
          & snowmelt(ic),                    & ! out
          & runoff_glac(ic),                 & ! out
          & pme_glacier = runoff(ic)         & ! out, P-E ... used as runoff for HD model
          & )
        infilt(ic)   = 0._wp
        drainage(ic) = 0._wp
        weq_snow_soil(ic) = 0._wp
        evapo_snow(ic) = evapopot(ic)
        IF (l_compat401_config) THEN
          snow_soil_dens(ic) = 330._wp
        ELSE
          snow_soil_dens(ic) = dens_snow
        END IF
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
    ELSE
      !$ACC DATA CREATE(skinres_max, skinres_canopy_max, weq_rootzone, weq_rootzone_max) &
      !$ACC   CREATE(trans_tmp, weq_soil_sl_tmp(1:nc,1:nsoil), weq_soil(1:nc)) &
      !$ACC   CREATE(weq_wsat_soil_sl_tmp(1:nc,1:nsoil), weq_wres_soil_sl_tmp(1:nc,1:nsoil)) &
      !$ACC   CREATE(ice_impedance_sl(1:nc,1:nsoil))

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO  ic = 1, nc
          weq_wsat_soil_sl_tmp(ic,is)      = vol_porosity_sl(ic,is)  * soil_depth_sl(ic,is)
          weq_wres_soil_sl_tmp(ic,is)      = vol_wres_sl(ic,is)      * soil_depth_sl(ic,is)
          ! @todo There should be no water in cells with field_cap == 0. Still, this is currently
          !       necessary. WHY?
          IF (weq_wsat_soil_sl_tmp(ic,is) > 0._wp) THEN
            weq_soil_sl_tmp(ic,is) = wtr_soil_sl(ic,is) + ice_soil_sl(ic,is)
          ELSE
            weq_soil_sl_tmp(ic,is) = 0._wp
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO  ic = 1, nc
        weq_soil(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP SEQ
      DO is = 1, nsoil
        !$ACC LOOP GANG VECTOR
        DO  ic = 1, nc
          weq_soil(ic) = weq_soil(ic) + weq_soil_sl_tmp(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL

      !$ACC WAIT(1)

      ! Compute actual and maximum total rootzone moisture as well as the total water amount in the soil (liquid + ice)
      ! Note, here weq_rootzone_max is computed based on porosity, not field capacity. This is because it is used to
      ! compute the bucket overflow in the ARNO scheme for infiltration and runoff, and therefore should consider the
      ! maximum amount of water that can fit into the rootzone.
      CALL get_amount_in_rootzone(weq_soil_sl_tmp(:,:),                                                  &
        &  dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl), dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl), &
        &  weq_rootzone(:))
      CALL get_amount_in_rootzone(weq_wsat_soil_sl_tmp(:,:),                                             &
        &  dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl), dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl), &
        &  weq_rootzone_max(:))

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        IF (tile%contains_vegetation) THEN
          skinres_canopy_max(ic) = w_skin_max_config * lai(ic) * fract_fpc_max(ic)
        ELSE
          skinres_canopy_max(ic) = 0._wp
        END IF
        IF (tile%contains_soil) THEN
          skinres_max(ic) =  w_skin_max_config + skinres_canopy_max(ic)
        ELSE
          skinres_max(ic) = 0._wp
        END IF
        IF (tile%is_vegetation) THEN
          trans_tmp(ic) = transpiration(ic)
        ELSE
          trans_tmp(ic) = 0._wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)

      ! Compute soil hydrological properties to get top layer ice impedance and
      ! soil hydrological conductivity for point scale infiltration
      CALL get_soilhyd_properties(                                                           &
        & dsl4jsb_Config(HYDRO_)%soilhydmodel,                                               &
        & dsl4jsb_Config(HYDRO_)%interpol_mean,                                              &
        & nc, nsoil, soil_depth_sl(:,:),                                                     &
        & wtr_soil_sl(:,:), ice_soil_sl(:,:), wtr_soil_pot_scool_sl(:,:),                &
        & weq_wsat_soil_sl_tmp(:,:), weq_wres_soil_sl_tmp(:,:),                              &
        & hyd_cond_sat_sl(:,:), matric_pot_sl(:,:), bclapp_sl(:,:), pore_size_index_sl(:,:), &
        & ice_impedance=ice_impedance_sl(:,:)                                                &
        )

      CALL calc_surface_hydrology_land(   &
        ! in
        & is_experiment_start,            & ! in
        & dtime,                          & ! in
        & ltpe_closed,                    & ! in
        & hydro_scale,                    & ! choice of hydrological scale scheme
        & l_dynsnow_config,               & ! in
        & l_infil_subzero_config,         & ! in
        & snow_depth_max_config,          & ! in
        & steepness,                      & ! Shape parameter describing orography    []
        & t_soil_sl(:,1),                 & ! Temperature of the uppermost soil layer [K]
        & t_unfilt(:),                    & ! Surface tempemperature [K] (unfiltered)
        & wind_10m(:),                    & ! wind speed at 10m height [m/s]
        & t_air(:),                       & ! lowest layer atmosphere temperature [K]
        & skinres_canopy_max(:),          & ! Capacity of the canopy skin reservoir (also used to limit snow on canopy) [m]
        & skinres_max(:),                 & ! Total capacity of the skin reservoirs, i.e. soil and canopy [m]
        & weq_pond_max(:),                & ! Total capacity of pond reservoir (water + ice) [m]
        & fract_snow(:),                  & ! surface snow fraction []
        & fract_skin(:),                  & ! skin reservoir fraction []
        & fract_pond(:),                  & ! Actual pond fraction []
        & hcap_grnd(:),                   & ! heat capacity of the uppermost soil layer [J m-2 K-1]
        & evapotrans(:),                  & ! evapotranspiration (including sublimation) [kg m-2 s-1]
        & evapopot(:),                    & ! potential evaporation (if there was enough water/ice) [kg m-2 s-1]
        & trans_tmp(:),                   & ! transpiration [kg m-2 s-1]
        & rain(:),                        & ! liquid precipitation [kg m-2 s-1]
        & snow(:),                        & ! solid precipitation [kg m-2 s-1]
        & weq_rootzone(:),                & ! total rootzone water content [m]
        & weq_rootzone_max(:),            & ! maximum total rootzone water content [m]
        & weq_soil(:),                    & ! total column soil moisture [m]
        & hyd_cond_sat_sl(:,1),           & ! saturated hydraulic conductivity of the top soil layer [m s-1]
        & ice_impedance_sl(:,1),          & ! ice impedance factor for top soil layer [-]
        ! inout
        & wtr_skin(:),                    & ! water content of the skin reservoir (vegetation and bare soil) [m]
        & weq_snow_soil(:),               & ! snow depth at the ground [m water equivalent]
        & weq_snow_can(:),                & ! snow depth on the canopy [m water equivalent]
        & snow_soil_dens(:),              & ! snow density on soil at non-glacier points [m water equivalent]
        & wtr_pond(:),                    & ! pond reservoir water content [m]
        & ice_pond(:),                    & ! pond reservoir ice content [m]
        ! out
        & q_snocpymlt(:),                 & ! Heating due to snow melt on canopy [W m-2]
        & snow_accum(:),                  & ! snow budget change within time step [m water equivalent]
        & snowmelt(:),                    & ! snow/ice melt at land points (excluding canopy) [kg m-2 s-1]
        & pond_freeze(:),                 & ! Amount of water freezing within ponds [kg m-2 s-1]
        & pond_melt(:),                   & ! Amount of ice melting within ponds [kg m-2 s-1]
        & evapotrans_soil(:),             & ! evapotranspiration from soil w/o skin, snow and pond reservoirs [kg m-2 s-1]
        & evapo_skin(:),                  & ! evaporation from skin reservoir [kg m-2 s-1]
        & evapo_snow(:),                  & ! evaporation from snow [kg m-2 s-1]
        & evapo_pond(:),                  & ! evaporation from ponds [kg m-2 s-1]
        & wtr_pond_net_flx(:),            & ! diagnostic net water inflow into surface water ponds [kg m-2 s-1]
        & water_to_soil(:),               & ! Water available for infiltration into the soil [m water equivalent]
        & evapo_deficit(:),               & ! Evaporation deficit flux due to inconsistent treatment of snow evap. [m]
        & infilt(:),                      & ! Infiltration flux into the soil [m water equivalent]
        & runoff(:),                      & ! Surface runoff [m water equivalent]
        & runoff_horton(:)                & ! Horton component of surface runoff [m water equivalent]
        & )

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      ! Update weq_pond storage
      !$ACC LOOP GANG VECTOR
      DO ic=1,nc
        weq_pond(ic) = wtr_pond(ic) + ice_pond(ic)
      END DO
      !$ACC END LOOP

      IF (l_compat401_config) THEN
        !$ACC LOOP GANG VECTOR
        DO ic=1,nc
          snow_soil_dens(ic) = 330._wp
        END DO
        !$ACC END LOOP
      END IF
      !$ACC END PARALLEL

      !$ACC WAIT(1)
      !$ACC END DATA
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
    IF (tile%is_glacier) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic=1,nc
          weq_snow(ic) = 0._wp
      END DO
      !$ACC END LOOP
    ELSE IF (tile%is_bare) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic=1,nc
        weq_snow(ic) = weq_snow_soil(ic)
      END DO
      !$ACC END LOOP
    ELSE IF (tile%contains_vegetation) THEN
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO ic=1,nc
        weq_snow(ic) = weq_snow_soil(ic) + weq_snow_can(ic)
      END DO
      !$ACC END LOOP
    END IF
    !$ACC END PARALLEL

    !$ACC WAIT(1)

  END SUBROUTINE update_surface_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_surface_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options


    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_hydrology'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_skin,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_soil,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_snow_can,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snowmelt,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, q_snocpymlt,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_glac,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_accum,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, snow_soil_dens,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans_soil,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_skin,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_snow,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_pond,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_to_soil,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, infilt,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff_horton,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_pond,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pond_freeze,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pond_melt,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond_net_flx, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_hydrology
  !
  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_soil_properties(tile, options)

    USE mo_hydro_constants, ONLY: vol_porosity_org_top, hyd_cond_sat_org_top, bclapp_org_top, matric_pot_org_top, &
      & pore_size_index_org_top, vol_field_cap_org_top, vol_p_wilt_org_top, vol_wres_org_top, &
      & vol_porosity_org_below, hyd_cond_sat_org_below, bclapp_org_below, matric_pot_org_below, &
      & pore_size_index_org_below, vol_field_cap_org_below, vol_p_wilt_org_below, vol_wres_org_below, &
      & thresh_org, beta_perc


    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & hyd_cond_sat,         &
      & vol_porosity,         &
      & bclapp,               &
      & matric_pot,           &
      & pore_size_index,      &
      & vol_field_cap,        &
      & vol_p_wilt,           &
      & vol_wres
    dsl4jsb_Real3D_onChunk :: &
      & fract_org_sl,         &
      & hyd_cond_sat_sl,      &
      & vol_porosity_sl,      &
      & bclapp_sl,            &
      & matric_pot_sl,        &
      & pore_size_index_sl,   &
      & vol_field_cap_sl,     &
      & vol_p_wilt_sl,        &
      & vol_wres_sl

    ! Locally allocated vectors
    !
    TYPE(t_jsb_model), POINTER :: model

    INTEGER  :: iblk, ics, ice, nc, ic
    INTEGER  :: nsoil, is

    REAL(wp) :: N_perc
    REAL(wp), ALLOCATABLE ::  fract_perc(:,:)          ! fraction of the soil with connected organic pathways
    REAL(wp), ALLOCATABLE ::  fract_uncon(:,:)         ! fraction of the soil with unconnected organic matter
    REAL(wp), ALLOCATABLE ::  hyd_cond_sat_uncon(:,:)  ! saturated hydraulic conductivity of 'unconnected' fraction

    REAL(wp), ALLOCATABLE ::  hyd_cond_sat_org(:,:)
    REAL(wp), ALLOCATABLE ::  vol_porosity_org(:,:)
    REAL(wp), ALLOCATABLE ::  bclapp_org(:,:)
    REAL(wp), ALLOCATABLE ::  matric_pot_org(:,:)
    REAL(wp), ALLOCATABLE ::  pore_size_index_org(:,:)
    REAL(wp), ALLOCATABLE ::  vol_field_cap_org(:,:)
    REAL(wp), ALLOCATABLE ::  vol_p_wilt_org(:,:)
    REAL(wp), ALLOCATABLE ::  vol_wres_org(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_properties'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%is_lake .OR. tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    nsoil = options%nsoil_w

    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(HYDRO_, hyd_cond_sat)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_porosity)       ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, bclapp)             ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, matric_pot)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, pore_size_index)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_field_cap)      ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_p_wilt)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, vol_wres)           ! in
    IF (dsl4jsb_Config(HYDRO_)%l_organic)   &
      &  dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)  ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_, hyd_cond_sat_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)    ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)          ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, matric_pot_sl)      ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, pore_size_index_sl) ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_cap_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_p_wilt_sl)      ! out
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_wres_sl)        ! out

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is=1,nsoil
      DO ic=1,nc
        hyd_cond_sat_sl   (ic,is) = hyd_cond_sat(ic)
        vol_porosity_sl   (ic,is) = vol_porosity(ic)
        bclapp_sl         (ic,is) = bclapp(ic)
        matric_pot_sl     (ic,is) = matric_pot(ic)
        pore_size_index_sl(ic,is) = pore_size_index(ic)
        vol_field_cap_sl  (ic,is) = vol_field_cap(ic)
        vol_p_wilt_sl     (ic,is) = vol_p_wilt(ic)
        vol_wres_sl       (ic,is) = vol_wres(ic)
      END DO
    END DO
    !$ACC END PARALLEL LOOP

    IF (dsl4jsb_Config(HYDRO_)%l_organic) THEN

      ALLOCATE(hyd_cond_sat_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(vol_porosity_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(bclapp_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(matric_pot_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(pore_size_index_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(vol_field_cap_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(vol_p_wilt_org(SIZE(hyd_cond_sat),nsoil))
      ALLOCATE(vol_wres_org(SIZE(hyd_cond_sat),nsoil))

      ALLOCATE(fract_perc(SIZE(hyd_cond_sat),nsoil))          ! fraction of the soil with connected organic pathways
      ALLOCATE(fract_uncon(SIZE(hyd_cond_sat),nsoil))         ! fraction of the soil with unconnected organic matter
      ALLOCATE(hyd_cond_sat_uncon(SIZE(hyd_cond_sat),nsoil))  ! saturated hydraulic conductivity of 'unconnected' fraction

      !$ACC DATA &
      !$ACC   CREATE(fract_perc, fract_uncon, hyd_cond_sat_uncon, hyd_cond_sat_org, vol_porosity_org) &
      !$ACC   CREATE(bclapp_org, matric_pot_org, pore_size_index_org, vol_field_cap_org, vol_p_wilt_org, vol_wres_org)

      ! Use different parameters for top layer and deeper layers
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        hyd_cond_sat_org(ic,1) = hyd_cond_sat_org_top
        vol_porosity_org(ic,1) = vol_porosity_org_top
        bclapp_org(ic,1) = bclapp_org_top
        matric_pot_org(ic,1) = matric_pot_org_top
        pore_size_index_org(ic,1) = pore_size_index_org_top
        vol_field_cap_org(ic,1) = vol_field_cap_org_top
        vol_p_wilt_org(ic,1) = vol_p_wilt_org_top
        vol_wres_org(ic,1) = vol_wres_org_top
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is=2,nsoil
        DO ic=1,nc
          hyd_cond_sat_org(ic,is) = hyd_cond_sat_org_below
          vol_porosity_org(ic,is) = vol_porosity_org_below
          bclapp_org(ic,is) = bclapp_org_below
          matric_pot_org(ic,is) = matric_pot_org_below
          pore_size_index_org(ic,is) = pore_size_index_org_below
          vol_field_cap_org(ic,is) = vol_field_cap_org_below
          vol_p_wilt_org(ic,is) = vol_p_wilt_org_below
          vol_wres_org(ic,is) = vol_wres_org_below
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ! Saturated hydraulic conductivity following percolation theory
      ! see CLM45 Tech Note, section 7.4.1 Hydraulic Properties, p.161 f.

      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is=1,nsoil
        DO ic=1,nc
          IF (fract_org_sl(ic,is) > thresh_org) THEN
            ! Connected flow pathways consisting of organic material only exist
            N_perc = (1._wp - thresh_org)**(-beta_perc)
            fract_perc(ic,is)  = MIN(1._wp, N_perc * fract_org_sl(ic,is) * (fract_org_sl(ic,is)-thresh_org)**beta_perc)
            fract_uncon(ic,is) = 1._wp - fract_perc(ic,is)
          ELSE
            ! No connected organic matter pathways exist, and flow passes mineral and organic soil components in series
            fract_perc(ic,is)  = 0._wp
            fract_uncon(ic,is) = 1._wp
          END IF

          IF (fract_org_sl(ic,is) > 0._wp) THEN
            hyd_cond_sat_uncon(ic,is) = fract_uncon(ic,is)                                           &
              &                       * ( (1._wp - fract_org_sl(ic,is)) / hyd_cond_sat_sl(ic,is)     &
              &                            + (fract_org_sl(ic,is) - fract_perc(ic,is)) / hyd_cond_sat_org(ic,is) )**(-1._wp)
            hyd_cond_sat_sl(ic,is)    = fract_uncon(ic,is) * hyd_cond_sat_uncon(ic,is)  &
              &                       + fract_perc(ic,is) * hyd_cond_sat_org(ic,is)
          END IF
        END DO
      END DO
      !$ACC END LOOP

      !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
      DO is=1,nsoil
        DO ic=1,nc
          vol_porosity_sl   (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_porosity_sl   (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_porosity_org(ic,is)
          bclapp_sl         (ic,is) = (1._wp - fract_org_sl(ic,is)) * bclapp_sl         (ic,is) &
            &                       + fract_org_sl(ic,is) * bclapp_org(ic,is)
          matric_pot_sl     (ic,is) = (1._wp - fract_org_sl(ic,is)) * matric_pot_sl     (ic,is) &
            &                       + fract_org_sl(ic,is) * matric_pot_org(ic,is)
          pore_size_index_sl(ic,is) = (1._wp - fract_org_sl(ic,is)) * pore_size_index_sl(ic,is) &
            &                       + fract_org_sl(ic,is) * pore_size_index_org(ic,is)
          vol_field_cap_sl  (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_field_cap_sl  (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_field_cap_org(ic,is)
          vol_p_wilt_sl     (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_p_wilt_sl     (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_p_wilt_org(ic,is)
          vol_wres_sl       (ic,is) = (1._wp - fract_org_sl(ic,is)) * vol_wres_sl       (ic,is) &
            &                       + fract_org_sl(ic,is) * vol_wres_org(ic,is)
        END DO
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

      !$ACC WAIT(1)
      !$ACC END DATA

      DEALLOCATE(hyd_cond_sat_org)
      DEALLOCATE(vol_porosity_org)
      DEALLOCATE(bclapp_org)
      DEALLOCATE(matric_pot_org)
      DEALLOCATE(pore_size_index_org)
      DEALLOCATE(vol_field_cap_org)
      DEALLOCATE(vol_p_wilt_org)
      DEALLOCATE(vol_wres_org)

      DEALLOCATE(fract_perc)
      DEALLOCATE(fract_uncon)
      DEALLOCATE(hyd_cond_sat_uncon)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_properties
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_properties(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_properties'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, hyd_cond_sat_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_porosity_sl,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, bclapp_sl,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, matric_pot_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, pore_size_index_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_field_cap_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_p_wilt_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, vol_wres_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_org_sl,       weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_properties
  !
  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "soil_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_soil_hydrology(tile, options)

    USE mo_hydro_process,          ONLY: calc_soil_hydrology

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: wtr_soilhyd_res
    dsl4jsb_Real2D_onChunk :: wtr_pond
    dsl4jsb_Real2D_onChunk :: ice_pond
    dsl4jsb_Real2D_onChunk :: weq_pond

    ! Locally allocated vectors

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid

    LOGICAL  :: ltpe_open
    LOGICAL  :: ltpe_closed

    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: dtime
    INTEGER  :: nsoil

    REAL(wp), POINTER :: area(:), lat(:), lon(:)
    REAL(wp) :: tile_fract(options%nc)
    LOGICAL :: l_fract(options%nc)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_hydrology'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    nsoil = options%nsoil_w
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%is_lake) RETURN

    ! Runoff and drainage for glaciers already computed in update_surface_hydrology
    ! Glaciers have no water in soil
    IF (tile%is_glacier) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    grid => Get_grid(model%grid_id)

    area => grid%area(ics:ice, iblk)
    lat => grid%lat(ics:ice, iblk)
    lon => grid%lon(ics:ice, iblk)

    dsl4jsb_Get_config(HYDRO_)

    ltpe_open = .FALSE.
    ltpe_closed = .FALSE.
    IF (model%config%tpe_scheme == 'open')   ltpe_open = .TRUE.
    IF (model%config%tpe_scheme == 'closed') ltpe_closed = .TRUE.

    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_soilhyd_res)  ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_, wtr_pond)         ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_, ice_pond)         ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_pond)         ! out

    !$ACC DATA CREATE(l_fract, tile_fract)

    CALL tile%Get_fraction(ics, ice, iblk, fract=tile_fract(:))
    l_fract(:) = tile_fract(:) > 0._wp
    !$ACC UPDATE DEVICE(l_fract, tile_fract) ASYNC(1)
    !$ACC WAIT(1)

    CALL calc_soil_hydrology( &
      ! in
      & nc, &
      & l_fract(:), &
      & lat(:), &
      & lon(:), &
      & nsoil, &
      & dtime, &
      & ltpe_closed, ltpe_open,                                 &
      & dsl4jsb_Config(HYDRO_)%enforce_water_budget,            & ! Water balance check setting
      & dsl4jsb_Config(HYDRO_)%soilhydmodel,                    &
      & dsl4jsb_Config(HYDRO_)%interpol_mean,                   &
      & dsl4jsb_Config(HYDRO_)%hydro_scale,                     & ! choice of hydrological scale scheme
      & dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract,               &
      & dsl4jsb_var3D_onChunk (HYDRO_, soil_depth_sl     ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, root_depth_sl     ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, hyd_cond_sat_sl   ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, matric_pot_sl     ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, bclapp_sl         ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, pore_size_index_sl),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, vol_porosity_sl   ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, vol_field_cap_sl  ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, vol_p_wilt_sl     ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, vol_wres_sl       ),     & ! volumetric residual soil water content [m m-1]
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pot_scool_sl ), &
      & dsl4jsb_var2D_onChunk (HYDRO_, fract_pond_max    ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, evapotrans_soil   ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, transpiration     ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, ice_pond          ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, weq_pond_max      ),     &
      ! inout
      & dsl4jsb_var2D_onChunk (HYDRO_, infilt            ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, runoff            ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_sl       ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, ice_soil_sl       ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond          ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, wtr_pond_net_flx  ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, tpe_overflow      ),     &
      & dsl4jsb_var2D_onChunk (HYDRO_, evapo_deficit     ),     &
      ! out
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_sat_sl   ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_fc_sl    ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_pwp_sl   ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_soil_res_sl   ),     & ! residual soil water content
      & dsl4jsb_var2D_onChunk (HYDRO_, runoff_dunne      ),     & ! dunne component of surface runoff
      & dsl4jsb_var2D_onChunk (HYDRO_, drainage          ),     &
      & dsl4jsb_var3D_onChunk (HYDRO_, drainage_sl       ),     & ! subsurface drainage on each soil layer
      & dsl4jsb_var3D_onChunk (HYDRO_, wtr_transp_down   ),     & ! lateral downards transport of water
      & dsl4jsb_var2D_onChunk (HYDRO_, wtr_soilhyd_res   )      &
      & )

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)

    ! Update weq_pond storage due to potential soil infiltration overflow
    !$ACC LOOP GANG VECTOR
    DO ic = 1, nc
      weq_pond(ic) = wtr_pond(ic) + ice_pond(ic)
    END DO

    ! Convert soil hydrology balance error to volume
    !$ACC LOOP GANG VECTOR
    DO ic = 1, nc
      wtr_soilhyd_res(ic) = wtr_soilhyd_res(ic) * tile_fract(ic) * area(ic)
    END DO
    !ACC END LOOP

    !$ACC END PARALLEL

    !$ACC WAIT(1)
    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_hydrology'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_soil_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_sat_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_fc_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_pwp_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_res_sl,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_pond,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_pond_net_flx, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, infilt,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, runoff_dunne,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, drainage_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_deficit,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soilhyd_res,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_transp_down,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_hydrology
  !
  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "canopy_cond_unstressed"
  !! Task "update_canopy_cond_unstressed" calculates canopy conductance without water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_canopy_cond_unstressed(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process, ONLY: get_canopy_conductance

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)
    dsl4jsb_Def_memory(A2L_)

    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: swpar_srf_down
    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited

    INTEGER :: iblk, ics, ice, nc, ic

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_unstressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(PHENO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      swpar_srf_down)        ! IN
    dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)                   ! IN
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_unlimited) ! OUT

    ! Compute (max) canopy (stomatal) conductance with no water stress
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      canopy_cond_unlimited(ic) = get_canopy_conductance( lai(ic), swpar_srf_down(ic) )
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE update_canopy_cond_unstressed
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "canopy_cond_unstressed"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_canopy_cond_unstressed(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_unstressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, canopy_cond_unlimited, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_unstressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "update" for task "water_stress"
  !! Task "update_water_stress" calculates water stress
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_water_stress(tile, options)
    ! ---------------------------
    ! Variables

    ! Used variables
    USE mo_hydro_process,      ONLY: get_water_stress_factor
    USE mo_hydro_util,         ONLY: get_amount_in_rootzone

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_config(SSE_)

    dsl4jsb_Def_memory(HYDRO_)

    ! HYDRO_
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: root_depth_sl
    dsl4jsb_Real3D_onChunk :: vol_field_cap_sl       !< Volumetric soil field capacity
    dsl4jsb_Real2D_onChunk :: wtr_rootzone           !< Liquid water content in the root zone
    dsl4jsb_Real2D_onChunk :: ice_rootzone           !< Frozen water content in the root zone
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl            !< Liquid water content in soil column
    dsl4jsb_Real3D_onChunk :: ice_soil_sl            !< Frozen water content in soil column
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl  !< Potentially suppercooled water content in soil column
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_scool_pot !< Potentially suppercooled water content in root zone
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_scool_act !< Actual supercooled water content in root zone
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_avail     !< Plant available liquid water content in root zone
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_avail_max !< Maximum plant available liquid water content in root zone
    dsl4jsb_Real2D_onChunk :: weq_rootzone_max       !< Maximum possible amount of water/ice in the rootzone
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: wtr_rootzone_rel       !< Relative root zone soil moisture

    INTEGER :: iblk, ics, ice, nc, ic, nsoil, is
    REAL(wp) :: config_w_soil_crit_fract, config_w_soil_wilt_fract
    REAL(wp) :: wtr_soil_pot_scool_sl_min(options%nc, options%nsoil_w)
    REAL(wp) :: weq_fcap_soil_sl(options%nc, options%nsoil_w)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_stress'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc
    nsoil = options%nsoil_w

    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_config(SSE_)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(HYDRO_)
    !
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    soil_depth_sl)          ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    vol_field_cap_sl)       ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    ice_soil_sl)            ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_pot_scool_sl)  ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone)           ! in
    dsl4jsb_Get_var3D_onChunk(HYDRO_,    root_depth_sl)          ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_rootzone_max)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)           ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    ice_rootzone)           ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_scool_pot) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_scool_act) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_avail)     ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_avail_max) ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_rootzone_rel)       ! out

    ! Get water stress

    ! Calculate water content and field capacity in the root zone
    ! @todo - do we need these?

    config_w_soil_crit_fract = dsl4jsb_Config(HYDRO_)%w_soil_crit_fract
    config_w_soil_wilt_fract = dsl4jsb_Config(HYDRO_)%w_soil_wilt_fract

    !$ACC DATA CREATE(wtr_soil_pot_scool_sl_min, weq_fcap_soil_sl)

    ! Calculate actual plant available water
    !   Meant here is all liquid water (without supercooled water) in the root zone. Water
    !   below the wilting point is not subtracted, so calling it available water might be misleading.

    ! Formerly, the root zone moisture was computed as part of the soil hydrology prior to diffusion
    !   and percolation. This seems wrong as the actual water stress should be the result of the
    !   updated soil state and not of the one at the start of the time step.
    CALL get_amount_in_rootzone(wtr_soil_sl(:,:),  &
                             &  soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone(:))

    ! Compute maximum root zone water content based on soil field capacity
    ! (replacing the one used for the ARNO scheme (surface hydrology) based on porosity)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
    DO is = 1, nsoil
      DO  ic = 1, nc
        weq_fcap_soil_sl(ic,is) = vol_field_cap_sl(ic,is) * soil_depth_sl(ic,is)
      END DO
    END DO
    !$ACC END PARALLEL LOOP
    CALL get_amount_in_rootzone(weq_fcap_soil_sl(:,:),  &
      &  soil_depth_sl(:,:), root_depth_sl(:,:), weq_rootzone_max(:))

    IF (dsl4jsb_Config(SSE_)%l_freeze .AND. dsl4jsb_Config(SSE_)%l_supercool) THEN
      ! Supercooled water is not available to plants
      CALL get_amount_in_rootzone(wtr_soil_pot_scool_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone_scool_pot(:))

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is=1,nsoil
        DO ic=1,nc
          wtr_soil_pot_scool_sl_min(ic,is) = MIN(wtr_soil_pot_scool_sl(ic,is), wtr_soil_sl(ic,is))
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      CALL get_amount_in_rootzone(wtr_soil_pot_scool_sl_min(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), wtr_rootzone_scool_act(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        wtr_rootzone_avail(ic) = MAX(wtr_rootzone(ic) - wtr_rootzone_scool_act(ic), 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      ! Available liquid water in the root zone
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        wtr_rootzone_avail(ic) = wtr_rootzone(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Calculate maximum plant available water
    !   As in the the available water calculation above, water below the wilting point is not
    !   subtracted.
    ! pdv: Here also the question is how to treat the potential rootzone soil moisture
    !   Definitly ice should be removed - otherwise plants suffer constant waterstress in permafrost-regions:
    !   If supercooled water is asssumed to be unavailable for plants then this should also be accounted for
    !   in wtr_rootzone_avail_max

    IF (dsl4jsb_Config(SSE_)%l_freeze .AND. dsl4jsb_Config(SSE_)%l_supercool) THEN
      CALL get_amount_in_rootzone(ice_soil_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        wtr_rootzone_avail_max(ic) = MAX(weq_rootzone_max(ic) - ice_rootzone(ic) - wtr_rootzone_scool_pot(ic), 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE IF (dsl4jsb_Config(SSE_)%l_freeze .AND. .NOT. dsl4jsb_Config(SSE_)%l_supercool) THEN
      CALL get_amount_in_rootzone(ice_soil_sl(:,:),  &
                               &  soil_depth_sl(:,:), root_depth_sl(:,:), ice_rootzone(:))
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        wtr_rootzone_avail_max(ic) = MAX(weq_rootzone_max(ic) - ice_rootzone(ic), 0._wp)
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        wtr_rootzone_avail_max(ic) = weq_rootzone_max(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    ! Calculate water stress
    ! @todo: maybe use the distributed permanent wilting point/field cap from above and change get_water_stress_factor function
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      water_stress(ic) = &
        & get_water_stress_factor(wtr_rootzone_avail(ic), wtr_rootzone_avail_max(ic), &
            &                 config_w_soil_crit_fract, config_w_soil_wilt_fract)
    END DO
    !$ACC END PARALLEL LOOP

    ! Calculate relative soil moisture (used in LoGro-P Phenology)
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic= 1, nc
      IF (wtr_rootzone_avail_max(ic) > 0._wp) THEN
        wtr_rootzone_rel(ic) =  MIN(wtr_rootzone_avail(ic) / wtr_rootzone_avail_max(ic), 1._wp)
      ELSE
        ! weq_rootzone_max is zero for glacier tiles
        wtr_rootzone_rel(ic) = 0._wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)
    !$ACC END DATA

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_stress

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "water_stress"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_water_stress(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(HYDRO_)

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_stress'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! HYDRO_
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, weq_rootzone_max,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_rootzone,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_rel,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_avail,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_rootzone_avail_max, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, water_stress,           weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_stress

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

    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)

    dsl4jsb_Real2D_onChunk :: canopy_cond_unlimited
    dsl4jsb_Real2D_onChunk :: canopy_cond_limited
    dsl4jsb_Real2D_onChunk :: water_stress
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: q_air
    dsl4jsb_Real2D_onChunk :: press_srf

    INTEGER :: iblk, ics, ice, nc, ic
    LOGICAL :: q_air_gt_qsat_tmp
    LOGICAL :: use_tmx

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_canopy_cond_stressed'

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice
    nc   = options%nc

    IF (.NOT. tile%Is_process_calculated(HYDRO_) .OR. .NOT. tile%is_vegetation) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    use_tmx = model%config%use_tmx

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      q_air)               ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      press_srf)           ! in
    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                   ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    water_stress)        ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_unlimited) ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    canopy_cond_limited) ! out

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
    !$ACC   PRIVATE(q_air_gt_qsat_tmp)
    DO ic=1,nc
      q_air_gt_qsat_tmp = q_air(ic) > qsat_water(t(ic),press_srf(ic), use_convect_tables=.NOT. use_tmx)
      ! Compute (actual) canopy (stomatal) conductance under water stress.
      canopy_cond_limited(ic) = get_canopy_conductance(canopy_cond_unlimited(ic), & ! in, unstressed canopy conductance
                                                       water_stress(ic),          & ! in, water stress factor
                                                       q_air_gt_qsat_tmp          & ! in, atmosphere saturated?
                                                      )
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

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

    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_canopy_cond_stressed'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, canopy_cond_limited, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_canopy_cond_stressed

  ! -------------------------------------------------------------------------------------------------------
  !>
  !> Implementation of "update" for task "evaporation"
  !! Task "update_evaporation" calculates evapo(transpi)ration
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_evaporation(tile, options)

    USE mo_phy_schemes,            ONLY: q_effective, heat_transfer_coef

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & q_acoef, &
      & q_bcoef, &
      & q_acoef_wtr, &
      & q_bcoef_wtr, &
      & q_acoef_ice, &
      & q_bcoef_ice, &
      & drag_srf, &
      & drag_wtr, &
      & drag_ice, &
      & ch, &
      & qsat_star, &
      & qsat_lwtr, &
      & qsat_lice, &
      & fact_qsat_srf, &
      & fact_qsat_trans_srf, &
      & fact_q_air, &
      & fract_lice, &
      & evapotrans_lnd, &
      & evapo_wtr, &
      & evapo_ice, &
      & evapopot,       &
      & evapotrans, &
      & transpiration

    REAL(wp) ::        &
      & q_air        , &  !< Humidity at lowest atmospheric level [kg kg-1]
      & q_air_eff    , &  !< Effective humidity at lowest atmospheric level
      & qsat_srf_eff , &  !< Effective surface saturation specific humidity
      & heat_tcoef        !< Heat transfer coefficient (rho*C_h*|v|)

    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: steplen, alpha

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_evaporation'

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    steplen = options%steplen
    alpha   = options%alpha

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(TURB_)

    ! Pointers to variables in memory
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapopot)   ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans) ! out

    IF (tile%contains_lake) THEN

      dsl4jsb_Get_var2D_onChunk(SEB_,     qsat_lwtr)    ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,   evapo_wtr)    ! out
      IF (model%config%use_tmx) THEN
        dsl4jsb_Get_var2D_onChunk(TURB_, ch)
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_acoef)       ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_bcoef)       ! in
      ELSE
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_acoef_wtr)   ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,  q_bcoef_wtr)   ! in
        dsl4jsb_Get_var2D_onChunk(A2L_,  drag_wtr)      ! in
      END IF

      IF (model%config%use_tmx) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
        !$ACC   PRIVATE(q_air, heat_tcoef)
        DO ic=1,nc
          q_air = q_bcoef(ic)                                   ! Old moisture at lowest atmospheric level
          heat_tcoef = ch(ic)                                   ! Transfer coefficient; TODO: distinguish wtr and ice?
          evapo_wtr (ic) = heat_tcoef * (q_air - qsat_lwtr(ic)) ! Potential evaporation
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
        !$ACC   PRIVATE(q_air, heat_tcoef)
        DO ic=1,nc
          q_air = q_acoef_wtr(ic) * qsat_lwtr(ic) + q_bcoef_wtr(ic) ! New moisture at lowest atmospheric level by back-substitution
          heat_tcoef = heat_transfer_coef(drag_wtr(ic), steplen, alpha)  ! Transfer coefficient
          evapo_wtr (ic) = heat_tcoef * (q_air - qsat_lwtr(ic))          ! Potential evaporation
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        IF (.NOT. model%config%use_tmx) THEN
          dsl4jsb_Get_var2D_onChunk(A2L_,   q_acoef_ice)       ! in
          dsl4jsb_Get_var2D_onChunk(A2L_,   q_bcoef_ice)       ! in
          dsl4jsb_Get_var2D_onChunk(A2L_,   drag_ice)    ! in
        END IF
        dsl4jsb_Get_var2D_onChunk(SEB_,   qsat_lice)         ! in
        dsl4jsb_Get_var2D_onChunk(SEB_,   fract_lice)        ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_, evapo_ice)         ! out

        IF (model%config%use_tmx) THEN
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
          !$ACC   PRIVATE(q_air, heat_tcoef)
          DO ic=1,nc
            q_air = q_bcoef(ic)                                   ! Old moisture at lowest atmospheric level
            heat_tcoef = ch(ic)                                   ! Transfer coefficient; TODO: distinguish between wtr and ice?
            evapo_ice (ic) = heat_tcoef * (q_air - qsat_lice(ic)) ! Potential evaporation
            evapopot  (ic) = (1._wp - fract_lice(ic)) * evapo_wtr(ic) + fract_lice(ic) * evapo_ice(ic)
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
          !$ACC   PRIVATE(q_air, heat_tcoef)
          DO ic=1,nc
            q_air = q_acoef_ice(ic) * qsat_lice(ic) + q_bcoef_ice(ic) ! New moisture at lowest atmospheric level by back-substitution
            heat_tcoef = heat_transfer_coef(drag_ice(ic), steplen, alpha)  ! Transfer coefficient
            evapo_ice (ic) = heat_tcoef * (q_air - qsat_lice(ic))          ! Potential evaporation
            evapopot  (ic) = (1._wp - fract_lice(ic)) * evapo_wtr(ic) + fract_lice(ic) * evapo_ice(ic)
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic=1,nc
          evapopot(ic) = evapo_wtr(ic)
        END DO
        !$ACC END PARALLEL LOOP
      END IF
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        evapotrans(ic) = evapopot(ic)
      END DO
      !$ACC END PARALLEL LOOP

    ELSE IF (tile%contains_soil .OR. tile%contains_glacier) THEN

      dsl4jsb_Get_var2D_onChunk(A2L_,  q_acoef)             ! in
      dsl4jsb_Get_var2D_onChunk(A2L_,  q_bcoef)             ! in
      dsl4jsb_Get_var2D_onChunk(SEB_,  qsat_star)           ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_q_air)          ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_srf)       ! in
      dsl4jsb_Get_var2D_onChunk(TURB_, fact_qsat_trans_srf) ! in

      IF (model%config%use_tmx) THEN
        dsl4jsb_Get_var2D_onChunk(TURB_, ch)
      ELSE
        dsl4jsb_Get_var2D_onChunk(A2L_,  drag_srf)          ! in
      END IF

      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans_lnd)           ! out
      IF (tile%contains_vegetation) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    transpiration)          ! out
      END IF

      IF (model%config%use_tmx) THEN
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
        !$ACC   PRIVATE(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
        DO ic=1,nc
          q_air = q_bcoef(ic)  ! Old moisture at lowest atmospheric level
          heat_tcoef = ch(ic)  ! Transfer coefficient

          ! Compute effective air moisture and surface saturation humidity
          q_air_eff      = q_effective( 0._wp, q_air, 1._wp, 0._wp)
          qsat_srf_eff   = q_effective(qsat_star(ic), q_air, fact_qsat_srf(ic), fact_q_air(ic))
          evapotrans_lnd(ic) = heat_tcoef * (q_air_eff - qsat_srf_eff)  ! Evapotranspiration
          evapotrans(ic)     = evapotrans_lnd(ic)
          evapopot(ic)       = heat_tcoef * (q_air     - qsat_star(ic)) ! Potential evaporation
          IF (tile%contains_vegetation) THEN
            transpiration(ic) = fact_qsat_trans_srf(ic) * evapopot(ic)  ! Transpiration
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
        !$ACC   PRIVATE(q_air, q_air_eff, qsat_srf_eff, heat_tcoef)
        DO ic=1,nc
          q_air = q_acoef(ic) * qsat_star(ic) + q_bcoef(ic)       ! New moisture at lowest atmospheric level by back-substitution
          heat_tcoef = heat_transfer_coef(drag_srf(ic), steplen, alpha)  ! Transfer coefficient

          ! Compute effective air moisture and surface saturation humidity
          q_air_eff      = q_effective( 0._wp, q_air, 1._wp, 0._wp)
          qsat_srf_eff   = q_effective(qsat_star(ic), q_air, fact_qsat_srf(ic), fact_q_air(ic))
          evapotrans_lnd(ic) = heat_tcoef * (q_air_eff - qsat_srf_eff)  ! Evapotranspiration
          evapotrans(ic)     = evapotrans_lnd(ic)
          evapopot(ic)       = heat_tcoef * (q_air     - qsat_star(ic)) ! Potential evaporation
          IF (tile%contains_vegetation) THEN
            transpiration(ic) = fact_qsat_trans_srf(ic) * evapopot(ic)  ! Transpiration
          END IF
        END DO
        !$ACC END PARALLEL LOOP
      END IF

    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

    !$ACC WAIT(1)

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_evaporation
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "evaporation"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_evaporation(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_evaporation'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapopot,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, transpiration,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapotrans_lnd,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_wtr,          weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(HYDRO_, evapo_ice,        weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_evaporation
  !
  !>
  !> Implementation of "update" for task "snow_and_ice_hydrology"
  !! Task "update_snow_hydrology" calculates the new snow depth and fraction
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_ice_hydrology'

    IF (options%nc > 0) CONTINUE ! avoid compiler warnings about dummy arguments not being used

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (tile%contains_lake) THEN
      !! Currently done in SEB process
    ELSE IF (tile%contains_land) THEN
!!$      CALL update_snow_and_ice_hydrology_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_snow_and_ice_hydrology
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snow_and_ice_hydrology"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_snow_and_ice_hydrology(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_ice_hydrology'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')


    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_ice_hydrology
  !
  ! ================================================================================================================================
  !
  !> Implementation of "update" for task "snow_and_wet_fraction"
  !! Task "update_snow_and_wet_fraction" calculates the new snow and wet fractions (skin and ponds)
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snow_and_wet_fraction( tile, options)

    ! Use declarations
    USE mo_hydro_process,          ONLY: calc_wskin_fractions_lice, calc_wet_fractions_veg, calc_wet_fractions_bare
    ! sollte jetzt in dieses File hier kommen: USE mo_hydro_process,  ONLY: calc_hydro_snow_and_skin_fraction
    USE mo_phy_schemes,            ONLY: heat_transfer_coef

    ! Arguments
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !----------------------------
    ! Local variables
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: grid
    REAL(wp), DIMENSION(options%nc) :: skinres_canopy_max, skinres_max, heat_tcoef

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snow_and_wet_fraction'

    INTEGER  :: iblk, ics, ice, nc, ic
    INTEGER  :: pond_dynamics_scheme
    REAL(wp) :: dtime, steplen, alpha, config_w_skin_max

    ! Declare pointers for process configuration and memory
    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(TURB_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(PHENO_)

    ! Declare pointers for variables in memory
    dsl4jsb_Real2D_onChunk :: lai
    dsl4jsb_Real2D_onChunk :: fract_snow
    dsl4jsb_Real2D_onChunk :: fract_skin
    dsl4jsb_Real2D_onChunk :: fract_wet
    dsl4jsb_Real2D_onChunk :: fract_pond
    dsl4jsb_Real2D_onChunk :: fract_pond_max
    dsl4jsb_Real2D_onChunk :: fract_fpc_max
    dsl4jsb_Real2D_onChunk :: fract_snow_can
    dsl4jsb_Real2D_onChunk :: fract_snow_soil
    dsl4jsb_Real2D_onChunk :: fract_snow_lice
    dsl4jsb_Real2D_onChunk :: weq_snow_can
    dsl4jsb_Real2D_onChunk :: weq_snow_soil
    dsl4jsb_Real2D_onChunk :: weq_snow_lice
    dsl4jsb_Real2D_onChunk :: weq_pond
    dsl4jsb_Real2D_onChunk :: weq_pond_max
    dsl4jsb_Real2D_onChunk :: wtr_skin
    dsl4jsb_Real2D_onChunk :: oro_stddev
    dsl4jsb_Real2D_onChunk :: drag_srf
    dsl4jsb_Real2D_onChunk :: ch
    dsl4jsb_Real2D_onChunk :: t
    dsl4jsb_Real2D_onChunk :: press_srf
    dsl4jsb_Real2D_onChunk :: q_air

    ! ---------------------------
    ! Go

    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    steplen = options%steplen
    alpha   = options%alpha

    IF (.NOT. tile%Is_process_calculated(HYDRO_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    ! Set pointers to process configs and memory
    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_config(HYDRO_)

    ! Use simple scalar for GPU
    config_w_skin_max = dsl4jsb_Config(HYDRO_)%w_skin_max
    pond_dynamics_scheme = dsl4jsb_Config(HYDRO_)%pond_dynamics

    dsl4jsb_Get_memory(HYDRO_)

    ! First handle LAKE_TYPE lct
    IF (tile%is_lake) THEN
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow_lice)    ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_lice)  ! out
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic=1,nc
          CALL calc_wskin_fractions_lice( &
            & weq_snow_lice(ic),          & ! in
            & fract_snow_lice(ic)         & ! out
            & )
        END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(1)
      END IF
      RETURN
    END IF

    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(SEB_)
    IF (tile%is_vegetation .OR. tile%is_land) THEN
      dsl4jsb_Get_memory(PHENO_)
    END IF

    ! Set pointers
    IF (model%config%use_tmx) THEN
      dsl4jsb_Get_memory(TURB_)
      dsl4jsb_Get_var2D_onChunk(TURB_, ch)            ! in
    ELSE
      dsl4jsb_Get_var2D_onChunk(A2L_, drag_srf)       ! in
    END IF
    dsl4jsb_Get_var2D_onChunk(A2L_, q_air)            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, press_srf)        ! in

    dsl4jsb_Get_var2D_onChunk(SEB_,      t)                ! in

    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    oro_stddev)        ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_pond_max)    ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_skin)          ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_pond_max)      ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_skin)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_pond)        ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_wet)         ! out
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_pond)          ! out
    END IF
    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_var2D_onChunk(PHENO_,    lai)               ! in
      dsl4jsb_Get_var2D_onChunk(PHENO_,    fract_fpc_max)     ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_can)    ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow_soil)    ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    fract_snow_soil)  ! out

#ifndef _OPENACC
    IF (ANY(t(:) < 50._wp .OR. t(:) > 400._wp)) THEN
      IF (ANY(t(:) < 50._wp)) THEN
        ic = MINLOC(t(:), DIM=1)
      ELSE
        ic = MaxLOC(t(:), DIM=1)
      END IF
      WRITE (message_text,*) 'Temperature out of bounds on tile ', tile%name, ' at ', '(', &
        &                    grid%lon(ic,iblk), ';', grid%lat(ic,iblk), '): t: ', t(ic)
      CALL finish(TRIM(routine), message_text)
    END IF
#endif

    !$ACC DATA CREATE(skinres_canopy_max, skinres_max, heat_tcoef)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic=1,nc
      ! Maximum capacity of skin reservoir (soil + canopy)
      IF (tile%contains_vegetation) THEN
        skinres_canopy_max(ic) = config_w_skin_max * lai(ic) * fract_fpc_max(ic)
      ELSE
        skinres_canopy_max(ic) = 0._wp
      END IF
      IF (tile%contains_soil) THEN
        skinres_max(ic) = config_w_skin_max + skinres_canopy_max(ic)
      ELSE
        skinres_max(ic) = 0._wp
      END IF

      IF (tile%is_glacier) THEN
        fract_snow(ic)      = 1._wp
        fract_snow_soil(ic) = 1._wp
      END IF
    END DO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)

    ! Transfer coefficient
    IF (model%config%use_tmx) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        heat_tcoef(ic) = ch(ic)  ! TODO: distinguish between wtr and ice?
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        heat_tcoef(ic) = heat_transfer_coef(drag_srf(ic), steplen, alpha)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (tile%contains_vegetation) THEN
      dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_can)       ! in
      CALL calc_wet_fractions_veg( &
        & dtime,                     & ! in
        & model%config%use_tmx,      & ! in
        & skinres_max(:),            & ! in
        & weq_pond_max(:),           & ! in
        & fract_pond_max(:),         & ! in
        & pond_dynamics_scheme,      & ! in
        & oro_stddev(:),             & ! in
        & t(:),                      & ! in, from the previous time step as long as this is called before the asselin filter
        & press_srf(:),              & ! in
        & heat_tcoef(:),             & ! in
        & q_air(:),                  & ! in
        & wtr_skin(:),               & ! in
        & weq_pond(:),               & ! in
        & weq_snow_soil(:),          & ! in
        & weq_snow_can(:),           & ! in
        & fract_snow_can(:),         & ! out
        & fract_skin(:),             & ! out
        & fract_pond(:),             & ! out
        & fract_wet(:),              & ! out
        & fract_snow_soil(:)         & ! out
        & )

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        fract_snow(ic) = fract_snow_soil(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (tile%is_bare) THEN
      CALL calc_wet_fractions_bare( &
        & dtime,                      & ! in
        & model%config%use_tmx,       & ! in
        & skinres_max,                & ! in
        & weq_pond_max(:),            & ! in
        & fract_pond_max(:),          & ! in
        & pond_dynamics_scheme,       & ! in
        & oro_stddev(:),              & ! in
        & t(:),                       & ! in
        & press_srf(:),               & ! in
        & heat_tcoef(:),              & ! in
        & q_air(:),                   & ! in
        & wtr_skin(:),                & ! in
        & weq_pond(:),                & ! in
        & weq_snow_soil(:),           & ! in
        & fract_skin(:),              & ! out
        & fract_pond(:),              & ! out
        & fract_wet(:),               & ! out
        & fract_snow_soil(:)          & ! out
        & )
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic=1,nc
        fract_snow(ic) = fract_snow_soil(ic)
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE update_snow_and_wet_fraction

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snow_and_wet_fraction"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_snow_and_wet_fraction(tile, options)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snow_and_wet_fraction'

    INTEGER  :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_wet,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_skin,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_pond,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow,       weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_soil,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_can,   weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(HYDRO_, fract_snow_lice,  weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snow_and_wet_fraction

  SUBROUTINE update_water_balance( tile, options)

    USE mo_jsb_physical_constants, ONLY: rhoh2o
    USE mo_jsb_impl_constants,     ONLY: WB_LOGGING, WB_ERROR

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER    :: model
    TYPE(t_jsb_grid),  POINTER    :: grid

    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)
    dsl4jsb_Def_memory(A2L_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & rain, &
      & snow, &
      & runoff, &
      & drainage, &
      & evapotrans, &
      & evapo_deficit, &
      & wtr_skin, &
      & weq_snow, &
      & weq_pond, &
      & weq_fluxes, &         ! [m3 s-1]
      & weq_land, &           ! [m3]
      & weq_balance_err, &    ! [m3/(time step)]
      & weq_balance_err_count ! [time steps]
    dsl4jsb_Real3D_onChunk :: &
      & wtr_soil_sl, &
      & ice_soil_sl

    REAL(wp), POINTER :: &
      & tile_fract(:), area(:), lat(:), lon(:)

    ! Locally allocated vectors
    !
    REAL(wp) :: &
      & tile_area,    &
      & wb_threshold, &
      & weq_land_new

    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: dtime
    LOGICAL  :: is_experiment_start, tile_contains_soil

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_water_balance'
    CHARACTER(len=4096)         :: message_text_long

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime
    is_experiment_start = is_time_experiment_start(options%current_datetime)

    ! IF (tile%lcts() RETURN   ! Check water balance only on root tile
    ! IF (tile%is_lake) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    grid => Get_grid(model%grid_id)

    tile_fract => tile%fract(ics:ice, iblk)
    area       => grid%area (ics:ice, iblk)
    lat        => grid%lat(ics:ice, iblk)
    lon        => grid%lon(ics:ice, iblk)

    tile_contains_soil = tile%contains_soil

    dsl4jsb_Get_config(HYDRO_)
    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_memory(A2L_)

    dsl4jsb_Get_var2D_onChunk(A2L_,      rain)                ! in
    dsl4jsb_Get_var2D_onChunk(A2L_,      snow)                ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    runoff)              ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    drainage)            ! in
    IF (tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_soil_sl)         ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_,    ice_soil_sl)         ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    wtr_skin)            ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_snow)            ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapo_deficit)       ! in
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_pond)            ! in
    END IF
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    evapotrans)            ! in
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_fluxes)            ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_land)              ! inout
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_balance_err)       ! out
    dsl4jsb_Get_var2D_onChunk(HYDRO_,    weq_balance_err_count) ! inout

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
    !$ACC   PRIVATE(tile_area, weq_land_new)
    DO ic=1,nc
      tile_area  = tile_fract(ic) * area(ic)

      weq_fluxes(ic) = rain(ic) + snow(ic) + evapotrans(ic) - runoff(ic) - drainage(ic)
      weq_fluxes(ic) = weq_fluxes(ic) * tile_area / rhoh2o               ! kg m-2 s-1 -> m3/s

      IF (tile_contains_soil) THEN
        weq_land_new = (wtr_skin(ic) + weq_snow(ic) + weq_pond(ic) + SUM(wtr_soil_sl(ic,:)) + SUM(ice_soil_sl(ic,:)))  &
          &              * tile_area          ! m -> m^3
      ELSE
        weq_land_new = 0._wp
      END IF

      IF (.NOT. is_experiment_start) THEN
        weq_balance_err(ic) = weq_land(ic) + weq_fluxes(ic) * dtime - weq_land_new
      END IF
      weq_land(ic) = weq_land_new

#ifndef _OPENACC
      ! Compute water balance threshold based on model resolution and time step. The threshold corresponds
      ! to a global water balance violation in the magnitude of 10cm of sea level rise after 1000 years
      wb_threshold = 1.0e-11_wp * tile_area * dtime

      IF (ABS(weq_balance_err(ic)) > wb_threshold) THEN
        weq_balance_err_count(ic) = weq_balance_err_count(ic) + 1.0_wp
        WRITE (message_text_long,*) 'Water balance violation [m3 dt-1]',               NEW_LINE('a'), &
          & 'on ',TRIM(tile%name),' tile at', lat(ic),'N and ',lon(ic),'E',            NEW_LINE('a'), &
          & '(ic: ',ic,' iblk: ',iblk, ' tile_fract:',tile_fract(ic),'):',             NEW_LINE('a'), &
          & 'WB Error:           ', weq_balance_err(ic),                               NEW_LINE('a'), &
          & 'Rainfall:           ', rain(ic)       * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Snowfall:           ', snow(ic)       * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Evapotranspiration: ', evapotrans(ic) * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Runoff:             ', runoff(ic)     * tile_area / rhoh2o * dtime,       NEW_LINE('a'), &
          & 'Drainage:           ', drainage(ic)   * tile_area / rhoh2o * dtime
        IF (tile%contains_soil) THEN
          WRITE (message_text_long,*) TRIM(message_text_long),                         NEW_LINE('a'), &
            & 'Skin reservoir:     ', wtr_skin(ic)           * tile_area,              NEW_LINE('a'), &
            & 'Snow reservoir:     ', weq_snow(ic)           * tile_area
          IF (.NOT. tile%is_lake) THEN
            WRITE (message_text_long,*) TRIM(message_text_long),                       NEW_LINE('a'), &
              & 'Pond reservoir:     ', weq_pond(ic)           * tile_area
          END IF
          WRITE (message_text_long,*) TRIM(message_text_long),                         NEW_LINE('a'), &
            & 'Soil water:         ', SUM(wtr_soil_sl(ic,:)) * tile_area,              NEW_LINE('a'), &
            & 'Soil ice:           ', SUM(ice_soil_sl(ic,:)) * tile_area,              NEW_LINE('a'), &
            & 'Old reservoirs:     ', weq_land(ic),                                    NEW_LINE('a'), &
            & 'ET deficit:         ', evapo_deficit(ic) * tile_area / rhoh2o * dtime
        END IF
        IF (dsl4jsb_Config(HYDRO_)%enforce_water_budget == WB_ERROR) THEN
          CALL finish (TRIM(routine), message_text_long)
        ELSE IF (dsl4jsb_Config(HYDRO_)%enforce_water_budget == WB_LOGGING) THEN
          CALL message (TRIM(routine), message_text_long, all_print=.TRUE.)
        END IF
      END IF
#endif
    END DO
    !$ACC END PARALLEL LOOP

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_water_balance

  SUBROUTINE aggregate_water_balance(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    !dsl4jsb_Def_memory(HYDRO_)

    !CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_water_balance'

    INTEGER :: iblk !, ics, ice

    iblk = options%iblk
    !ics  = options%ics
    !ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    !dsl4jsb_Get_memory(HYDRO_)

    ! weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    ! Don't aggregate, but explicitely compute water balance on each tile
    CALL update_water_balance(tile, options)
    ! dsl4jsb_Aggregate_onChunk(HYDRO_, weq_fluxes,       weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(HYDRO_, weq_land,         weighted_by_fract)
    ! dsl4jsb_Aggregate_onChunk(HYDRO_, weq_balance_err,  weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_water_balance

  !-----------------------------------------------------------------------------------------------------
  !> Global land mean hydrology output
  !!
  !! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_hydrology_diagnostics(tile)
#ifndef __ICON__
    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_hydrology_diagnostics'
    IF (debug_on()) CALL message(TRIM(routine), 'Global diagnostics only available with ICON')
#else

    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid
    USE mo_jsb_impl_constants,    ONLY: WB_IGNORE

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_hydrology_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: transpiration
    dsl4jsb_Real2D_onDomain :: evapotrans
    dsl4jsb_Real2D_onDomain :: weq_land
    dsl4jsb_Real2D_onDomain :: discharge_ocean
    dsl4jsb_Real2D_onDomain :: wtr_rootzone_rel
    dsl4jsb_Real2D_onDomain :: fract_snow
    dsl4jsb_Real2D_onDomain :: weq_snow
    dsl4jsb_Real2D_onDomain :: weq_balance_err
    dsl4jsb_Real2D_onDomain :: weq_balance_err_count

    LOGICAL, SAVE           :: print_wb_warning = .TRUE.

    REAL(wp), POINTER       :: trans_gmean(:)
    REAL(wp), POINTER       :: evapotrans_gmean(:)
    REAL(wp), POINTER       :: weq_land_gsum(:)
    REAL(wp), POINTER       :: discharge_ocean_gsum(:)
    REAL(wp), POINTER       :: wtr_rootzone_rel_gmean(:)
    REAL(wp), POINTER       :: fract_snow_gsum(:)
    REAL(wp), POINTER       :: weq_snow_gsum(:)
    REAL(wp), POINTER       :: weq_balance_err_gsum(:)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)
    REAL(wp)               :: global_land_area


    dsl4jsb_Get_memory(HYDRO_)
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  transpiration)              ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  evapotrans)                 ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_land)                   ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  discharge_ocean)            ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  wtr_rootzone_rel)           ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  fract_snow)                 ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_snow)                   ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_balance_err)            ! in
    dsl4jsb_Get_var2D_onDomain(HYDRO_,  weq_balance_err_count)      ! in


    trans_gmean          => HYDRO__mem%trans_gmean%ptr(:)           ! out
    evapotrans_gmean     => HYDRO__mem%evapotrans_gmean%ptr(:)      ! out
    weq_land_gsum        => HYDRO__mem%weq_land_gsum%ptr(:)         ! out
    discharge_ocean_gsum => HYDRO__mem%discharge_ocean_gsum%ptr(:)  ! out
    wtr_rootzone_rel_gmean => HYDRO__mem%wtr_rootzone_rel_gmean%ptr(:)  ! out
    fract_snow_gsum      => HYDRO__mem%fract_snow_gsum%ptr(:)       ! out
    weq_snow_gsum        => HYDRO__mem%weq_snow_gsum%ptr(:)         ! out
    weq_balance_err_gsum => HYDRO__mem%weq_balance_err_gsum%ptr(:)  ! out


    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea

    dsl4jsb_Get_config(HYDRO_)

    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')


    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    ! Domain Mask - to mask all halo cells for global sums (otherwise these
    ! cells are counted twice)
    ALLOCATE (in_domain(grid%nproma,grid%nblks))
    WHERE (is_in_domain(:,:))
      in_domain = 1._wp
    ELSEWHERE
      in_domain = 0._wp
    END WHERE

    ALLOCATE (scaling(grid%nproma,grid%nblks))

    ! Calculate 1d global land variables, if requested for output
    global_land_area = global_sum_array(area(:,:) * notsea(:,:) * in_domain(:,:))
    scaling(:,:) = notsea(:,:) * area(:,:) * in_domain(:,:)
    IF (HYDRO__mem%trans_gmean%is_in_output)           &
      &  trans_gmean          = global_sum_array(transpiration(:,:)   * scaling(:,:)) / global_land_area
    IF (HYDRO__mem%evapotrans_gmean%is_in_output)      &
      &  evapotrans_gmean     = global_sum_array(evapotrans(:,:)      * scaling(:,:)) / global_land_area
    ! Unit transformation from [m3] to [km3]: 1.e-9
    IF (HYDRO__mem%weq_land_gsum%is_in_output)         &
      &  weq_land_gsum        = global_sum_array(weq_land(:,:)  * notsea(:,:) * in_domain(:,:)) * 1.e-9_wp
    ! Unit transformation from [m3/s] to [Sv] (1 Sv = 1.e6 m3/s)
    IF (HYDRO__mem%discharge_ocean_gsum%is_in_output)  &
      &  discharge_ocean_gsum = global_sum_array(discharge_ocean(:,:) * in_domain(:,:)) * 1.e-6_wp
    IF (HYDRO__mem%wtr_rootzone_rel_gmean%is_in_output)  &
      &  wtr_rootzone_rel_gmean = global_sum_array(wtr_rootzone_rel(:,:) * scaling(:,:)) / global_land_area
    ! Unit transformation from [m2] to [Mio km2]: 1.e-12
    IF (HYDRO__mem%fract_snow_gsum%is_in_output)       &
      &  fract_snow_gsum      = global_sum_array(fract_snow(:,:)      * scaling(:,:)) * 1.e-12_wp
    ! Unit transformation from [m water equivalent](= [t]) to [Gt]: 1.e-9
    IF (HYDRO__mem%weq_snow_gsum%is_in_output)         &
      &  weq_snow_gsum        = global_sum_array(weq_snow(:,:)        * scaling(:,:)) * 1.e-9_wp
    ! No unit transformation [m3]: 1
    IF (HYDRO__mem%weq_balance_err_gsum%is_in_output)  &
      &  weq_balance_err_gsum = global_sum_array(weq_balance_err(:,:) * notsea(:,:) * in_domain(:,:))

    IF (dsl4jsb_Config(HYDRO_)%enforce_water_budget == WB_IGNORE .AND. print_wb_warning) THEN
      IF (global_sum_array(weq_balance_err_count(:,:) * in_domain(:,:)) > 0.5_wp) THEN
        CALL message (TRIM(routine), 'Water balance issues detected. '// &
          & 'Consider rerun the simulation with enforce_water_budget = "logging" for more details.')
        ! Don't repeat this warning during this simulation period
        print_wb_warning = .FALSE.
      END IF
    END IF

    DEALLOCATE (scaling, in_domain)
#endif
  END SUBROUTINE global_hydrology_diagnostics

#endif
END MODULE mo_hydro_interface
