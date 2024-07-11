!> Contains the interfaces for the surface energy balance process.
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

MODULE mo_seb_interface
#ifndef __NO_JSBACH__

  ! -------------------------------------------------------------------------------------------------------
  ! Used variables of module

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_jsb_time,        ONLY: is_time_experiment_start
  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message, finish

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SEB_, HYDRO_, SSE_
  dsl4jsb_Use_config(SEB_)
  dsl4jsb_Use_config(SSE_)

  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)
  dsl4jsb_Use_memory(SSE_)

  ! -------------------------------------------------------------------------------------------------------
  ! Module variables

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_seb_tasks, global_seb_diagnostics

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_energy
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_energy
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_energy
  END TYPE tsk_surface_energy

  INTERFACE tsk_surface_energy
    PROCEDURE Create_task_surface_energy
  END INTERFACE tsk_surface_energy

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_snowmelt_correction
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_snowmelt_correction
    PROCEDURE, NOPASS :: Aggregate => aggregate_snowmelt_correction
  END TYPE tsk_snowmelt_correction

  INTERFACE tsk_snowmelt_correction
    PROCEDURE Create_task_snowmelt_correction
  END INTERFACE tsk_snowmelt_correction

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_asselin_filter
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_asselin
    PROCEDURE, NOPASS :: Aggregate => aggregate_asselin
  END TYPE tsk_asselin_filter

  INTERFACE tsk_asselin_filter
    PROCEDURE Create_task_asselin_filter
  END INTERFACE tsk_asselin_filter

  TYPE, EXTENDS(t_jsb_process_task) :: tsk_surface_fluxes
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_surface_fluxes
    PROCEDURE, NOPASS :: Aggregate => aggregate_surface_fluxes
  END TYPE tsk_surface_fluxes

  INTERFACE tsk_surface_fluxes
    PROCEDURE Create_task_surface_fluxes
  END INTERFACE tsk_surface_fluxes

  CHARACTER(len=*), PARAMETER :: modname = 'mo_seb_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_energy task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_energy"
  !!
  FUNCTION Create_task_surface_energy(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_energy::return_ptr)
    CALL return_ptr%Construct(name='surface_energy', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for snowmelt_correction task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "snowmelt_correction"
  !!
  FUNCTION Create_task_snowmelt_correction(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_snowmelt_correction::return_ptr)
    CALL return_ptr%Construct(name='snowmelt_correction', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_snowmelt_correction

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for asselin_filter task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "asselin_filter"
  !!
  FUNCTION Create_task_asselin_filter(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_asselin_filter::return_ptr)
    CALL return_ptr%Construct(name='asselin_filter', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_asselin_filter

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for surface_energy fluxes
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "surface_fluxes"
  !!
  FUNCTION Create_task_surface_fluxes(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_surface_fluxes::return_ptr)
    CALL return_ptr%Construct(name='surface_fluxes', process_id=SEB_, owner_model_id=model_id)

  END FUNCTION Create_task_surface_fluxes

  ! =======================================================================================================
  !> Register tasks for surface energy process
  !!
  !! @param[in,out] this      Instance of surface energy process class
  !! @param[in]     model_id  Model id
  !!
  SUBROUTINE Register_seb_tasks(this, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: this
    INTEGER,              INTENT(in)    :: model_id

    CALL this%Register_task(tsk_surface_energy(model_id))
    CALL this%Register_task(tsk_snowmelt_correction(model_id))
    CALL this%Register_task(tsk_asselin_filter(model_id))
    CALL this%Register_task(tsk_surface_fluxes(model_id))

  END SUBROUTINE Register_seb_tasks

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "surface energy"
  !! Task "update_surface_energy" calculates the new surface temperature from the surface energy balance (sensible and latent heat,
  !! net radiation, ground heat).
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_energy(tile, options)

    USE mo_seb_land, ONLY: update_surface_energy_land
    USE mo_seb_lake, ONLY: update_surface_energy_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_energy'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%is_lake) THEN
      CALL update_surface_energy_lake(tile, options)
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_surface_energy_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type '//TRIM(tile%lcts(1)%name)//' on tile '//TRIM(tile%name))
    END IF

  END SUBROUTINE update_surface_energy

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_energy"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_surface_energy(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_config(SEB_)

    TYPE(t_jsb_model), POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_energy'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SEB_, t,             weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_old,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_unfilt,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_eff4,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, qsat_star,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, s_star,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, forc_hflx,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, heat_cap,      weighted_by_fract)

    IF (model%config%use_lakes) THEN
      dsl4jsb_Aggregate_onChunk(SEB_, t_lwtr,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, s_lwtr,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, qsat_lwtr,  weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_, fract_lice, weighted_by_fract)
      IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
        dsl4jsb_Aggregate_onChunk(SEB_, t_lice,     weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(SEB_, s_lice,     weighted_by_fract)
        dsl4jsb_Aggregate_onChunk(SEB_, qsat_lice,  weighted_by_fract)
      END IF
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_energy

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "snowmelt_correction"
  !! Task "snowmelt_correction" applies a correction for t_unfilt due to snow and ice melting,
  !! surface water freezing and thawing and soil moisture freezing and thawing.
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_snowmelt_correction(tile, options)

    USE mo_jsb_physical_constants, ONLY: alf

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: snowmelt
    dsl4jsb_Real2D_onChunk :: pond_freeze
    dsl4jsb_Real2D_onChunk :: pond_melt
    dsl4jsb_Real2D_onChunk :: hcap_grnd_old
    dsl4jsb_Real3D_onChunk :: snow_depth_sl
    dsl4jsb_Real3D_onChunk :: vol_heat_cap_sl
    dsl4jsb_Real3D_onChunk :: wtr_freeze_sl
    dsl4jsb_Real3D_onChunk :: ice_melt_sl


    INTEGER  :: iblk , ics, ice, nc, ic
    REAL(wp) :: dtime, dz1

    LOGICAL :: &
      & l_freeze_config, &
      & l_snow_config
    REAL(wp), DIMENSION(options%nc) :: &
      & heat_cap

    TYPE(t_jsb_vgrid), POINTER :: soil_e
    INTEGER                    :: nsnow

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_snowmelt_correction'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    nsnow = options%nsnow_e
    dtime = options%dtime

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SSE_)
    l_freeze_config = dsl4jsb_Config(SSE_)%l_freeze
    l_snow_config   = dsl4jsb_Config(SSE_)%l_snow
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)
    IF (tile%contains_soil .OR. tile%is_glacier) THEN
      dsl4jsb_Get_memory(SSE_)
    END IF

    soil_e => Get_vgrid('soil_depth_energy')
    dz1 = soil_e%dz(1)

    IF (tile%contains_lake) THEN
      ! No correction for lake tile

    ELSE IF (tile%contains_soil .OR. tile%is_glacier) THEN

      dsl4jsb_Get_var2D_onChunk(SEB_,      t_unfilt)         ! inout
      dsl4jsb_Get_var2D_onChunk(HYDRO_,    snowmelt)         ! in
      dsl4jsb_Get_var2D_onChunk(SSE_,      hcap_grnd_old)    ! in
      IF (dsl4jsb_Config(SSE_)%l_snow) THEN
        dsl4jsb_Get_var3D_onChunk(SSE_,      snow_depth_sl)    ! in
      END IF
      IF (.NOT. tile%is_glacier) THEN
        dsl4jsb_Get_var3D_onChunk(SSE_,      vol_heat_cap_sl) ! in
        dsl4jsb_Get_var3D_onChunk(HYDRO_,    wtr_freeze_sl)   ! in
        dsl4jsb_Get_var3D_onChunk(HYDRO_,    ice_melt_sl)     ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    pond_freeze)     ! in
        dsl4jsb_Get_var2D_onChunk(HYDRO_,    pond_melt)       ! in
      END IF

      IF (.NOT. is_time_experiment_start(options%current_datetime)) THEN
        !$ACC DATA CREATE(heat_cap) ASYNC(1)
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic=1,nc

          ! Correction for snowmelt and pond melting/freezing
          t_unfilt(ic) = t_unfilt(ic) - snowmelt(ic)    * dtime * alf / hcap_grnd_old(ic)
          IF (.NOT. tile%is_glacier) THEN
            t_unfilt(ic) = t_unfilt(ic) - pond_melt(ic)   * dtime * alf / hcap_grnd_old(ic)
            t_unfilt(ic) = t_unfilt(ic) + pond_freeze(ic) * dtime * alf / hcap_grnd_old(ic)

            IF (l_freeze_config) THEN
              heat_cap(ic) = vol_heat_cap_sl(ic,1) * dz1
              IF (l_snow_config) THEN
                IF (snow_depth_sl(ic,nsnow) < EPSILON(1._wp) .AND. heat_cap(ic) > 0._wp) THEN  ! No snow layers present
                  ! Correction for melting of ice in soil layers
                  t_unfilt(ic) = t_unfilt(ic) - ice_melt_sl(ic,1)   * dtime * alf / heat_cap(ic)
                  t_unfilt(ic) = t_unfilt(ic) + wtr_freeze_sl(ic,1) * dtime * alf / heat_cap(ic)
                END IF
              ELSE   ! No snow layers present
                IF (heat_cap(ic) > 0._wp) THEN
                  ! Correction for melting of ice in soil layers
                  t_unfilt(ic) = t_unfilt(ic) - ice_melt_sl(ic,1)   * dtime * alf / heat_cap(ic)
                  t_unfilt(ic) = t_unfilt(ic) + wtr_freeze_sl(ic,1) * dtime * alf / heat_cap(ic)
                END IF
              END IF
            END IF
          END IF

        END DO
        !$ACC END PARALLEL LOOP
        !$ACC WAIT(1)
        !$ACC END DATA
      END IF

    ELSE

      CALL finish(TRIM(routine), 'Called for invalid lct_type')

    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_snowmelt_correction

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "snowmelt_correction"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_snowmelt_correction(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_snowmelt_correction'

    INTEGER :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SEB_, t_unfilt, weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_snowmelt_correction

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "asselin_filter"
  !! Task "asselin_filter" calculates applies the Asselin time filter to the new surface temperature, if applicable
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     config  Vector of process configurations.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_asselin(tile, options)

    USE mo_seb_land, ONLY: update_asselin_land

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_asselin'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%contains_lake) THEN
      ! No Asselin filter for lake tile
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_asselin_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_asselin_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_asselin

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "asselin_filter"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_asselin(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SEB_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_asselin'

    INTEGER :: iblk , ics, ice
    TYPE(t_jsb_model), POINTER :: model

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SEB_, t,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_, t_filt, weighted_by_fract)
    IF (model%config%use_tmx) THEN
      dsl4jsb_Aggregate_onChunk(SEB_, t_rad4, weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_asselin

  ! ================================================================================================================================
  !>
  !> Implementation of "update" for task "surface fluxes"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE update_surface_fluxes(tile, options)

    USE mo_seb_land, ONLY: update_surface_fluxes_land
    USE mo_seb_lake, ONLY: update_surface_fluxes_lake

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_surface_fluxes'

    IF (.NOT. tile%Is_process_calculated(SEB_)) RETURN

    model => Get_model(tile%owner_model_id)

    IF (tile%contains_lake) THEN
      CALL update_surface_fluxes_lake(tile, options)
    ELSE IF (tile%is_land .AND. .NOT. model%config%use_tmx) THEN
      CALL update_surface_fluxes_land(tile, options)
    ELSE IF (model%config%use_tmx) THEN
      CALL update_surface_fluxes_land(tile, options)
    ELSE
      CALL finish(TRIM(routine), 'Called for invalid lct_type')
    END IF

  END SUBROUTINE update_surface_fluxes

  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "surface_fluxes"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  SUBROUTINE aggregate_surface_fluxes(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_config(SEB_)
    dsl4jsb_Def_memory(SEB_)

    TYPE(t_jsb_model),       POINTER :: model
    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_surface_fluxes'

    INTEGER :: iblk , ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SEB_)
    dsl4jsb_Get_memory(SEB_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx,      weighted_by_fract)

    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_lnd,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_lnd,  weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_wtr,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_wtr,  weighted_by_fract)
    IF (dsl4jsb_Config(SEB_)%l_ice_on_lakes) THEN
      dsl4jsb_Aggregate_onChunk(SEB_,   latent_hflx_ice,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SEB_,   sensible_hflx_ice,  weighted_by_fract)
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_surface_fluxes

  !-----------------------------------------------------------------------------------------------------
  !> Calculations of diagnostic global land mean output
  !!
  !! The routine is called from jsbach_finish_timestep, after the loop over the nproma blocks.
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE global_seb_diagnostics(tile)

#ifndef __ICON__
    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_seb_diagnostics'
    IF (debug_on()) CALL message(TRIM(routine), 'Global diagnostics only available with ICON')
#else

    USE mo_sync,                  ONLY: global_sum_array
    USE mo_jsb_grid,              ONLY: Get_grid
    USE mo_jsb_grid_class,        ONLY: t_jsb_grid

    ! Argument
    CLASS(t_jsb_tile_abstract), INTENT(in) :: tile

    ! Local variables
    !
    dsl4jsb_Def_memory(SEB_)

    CHARACTER(len=*),  PARAMETER  :: routine = modname//':global_seb_diagnostics'

    ! Pointers to variables in memory

    dsl4jsb_Real2D_onDomain :: t

    REAL(wp), POINTER       :: t_gmean(:)

    TYPE(t_jsb_model), POINTER      :: model
    TYPE(t_jsb_grid),  POINTER      :: grid

    REAL(wp), POINTER      :: area(:,:)
    REAL(wp), POINTER      :: notsea(:,:)
    LOGICAL,  POINTER      :: is_in_domain(:,:) ! T: cell in domain (not halo)
    REAL(wp), ALLOCATABLE  :: in_domain (:,:)   ! 1: cell in domain, 0: halo cell
    REAL(wp), ALLOCATABLE  :: scaling (:,:)
    REAL(wp)               :: global_land_area
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_var2D_onDomain(SEB_,  t)        ! in

    t_gmean        => SEB__mem%t_gmean%ptr(:)   ! out

    model => Get_model(tile%owner_model_id)
    grid  => Get_grid(model%grid_id)
    area         => grid%area(:,:)
    is_in_domain => grid%patch%cells%decomp_info%owner_mask(:,:)
    notsea       => tile%fract(:,:)   ! fraction of the box tile: notsea


    IF (debug_on()) CALL message(TRIM(routine), 'Starting routine')


    IF (ASSOCIATED(tile%parent)) CALL finish(TRIM(routine), 'Should only be called for the root tile')

    IF (SEB__mem%t_gmean%is_in_output) THEN

      ! Domain Mask - to mask all halo cells for global sums (otherwise these
      ! cells are counted twice)
      ALLOCATE (in_domain(grid%nproma,grid%nblks))
      WHERE (is_in_domain(:,:))
        in_domain = 1._wp
      ELSEWHERE
        in_domain = 0._wp
      END WHERE

      ALLOCATE (scaling(grid%nproma,grid%nblks))

      ! Calculate global land mean seb variables
      global_land_area = global_sum_array(area(:,:) * notsea(:,:) * in_domain(:,:))
      scaling(:,:) = notsea(:,:) * area(:,:) * in_domain(:,:)
      t_gmean      = global_sum_array(t(:,:) * scaling(:,:)) / global_land_area

      DEALLOCATE (scaling, in_domain)
    END IF

#endif
  END SUBROUTINE global_seb_diagnostics

  !
  ! ================================================================================================================================
  !
#endif
END MODULE mo_seb_interface
