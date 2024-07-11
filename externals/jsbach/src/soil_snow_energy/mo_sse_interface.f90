!> Contains the interfaces to the soil and snow energy processes
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
!NEC$ options "-finline-file=externals/jsbach/src/soil_snow_energy/mo_sse_process.pp-jsb.f90"

MODULE mo_sse_interface
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message

  USE mo_jsb_control,     ONLY: debug_on
  USE mo_jsb_time,        ONLY: is_time_experiment_start

  USE mo_jsb_model_class,    ONLY: t_jsb_model
  USE mo_jsb_class,          ONLY: Get_model
  USE mo_jsb_grid_class,     ONLY: t_jsb_vgrid
  USE mo_jsb_grid,           ONLY: Get_vgrid
  USE mo_jsb_tile_class,     ONLY: t_jsb_tile_abstract, t_jsb_aggregator
  USE mo_jsb_process_class,  ONLY: t_jsb_process
  USE mo_jsb_task_class,     ONLY: t_jsb_process_task, t_jsb_task_options

  dsl4jsb_Use_processes SSE_, SEB_, HYDRO_
  dsl4jsb_Use_config(SSE_)
  dsl4jsb_Use_config(HYDRO_)
  dsl4jsb_Use_memory(SSE_)
  dsl4jsb_Use_memory(SEB_)
  dsl4jsb_Use_memory(HYDRO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Register_sse_tasks

  !> Type definition for soil_temp_and_hydro task
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_and_snow_temp
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_and_snow_temperature
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_and_snow_temperature
  END TYPE tsk_soil_and_snow_temp

  !> Constructor interface for soil_temperature task
  INTERFACE tsk_soil_and_snow_temp
    PROCEDURE Create_task_soil_and_snow_temperature
  END INTERFACE tsk_soil_and_snow_temp

  !> Type definition for soil_properties
  TYPE, EXTENDS(t_jsb_process_task) :: tsk_soil_and_snow_properties
  CONTAINS
    PROCEDURE, NOPASS :: Integrate => update_soil_and_snow_properties
    PROCEDURE, NOPASS :: Aggregate => aggregate_soil_and_snow_properties
  END TYPE tsk_soil_and_snow_properties

  !> Constructor interface for soil_properties
  INTERFACE tsk_soil_and_snow_properties
    PROCEDURE Create_task_soil_and_snow_properties
  END INTERFACE tsk_soil_and_snow_properties

  CHARACTER(len=*), PARAMETER :: modname = 'mo_sse_interface'

CONTAINS

  ! ================================================================================================================================
  !! Constructors for tasks

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_temperature
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_temperature"
  !!
  FUNCTION Create_task_soil_and_snow_temperature(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_and_snow_temp::return_ptr)
    CALL return_ptr%Construct(name='soil_and_snow_temperature', process_id=SSE_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_and_snow_temperature

  ! -------------------------------------------------------------------------------------------------------
  !> Constructor for soil_properties task
  !!
  !! @param[in]     model_id     Model id
  !! @return        return_ptr   Instance of process task "soil_properties"
  !!
  FUNCTION Create_task_soil_and_snow_properties(model_id) RESULT(return_ptr)

    INTEGER,                   INTENT(in) :: model_id
    CLASS(t_jsb_process_task), POINTER    :: return_ptr

    ALLOCATE(tsk_soil_and_snow_properties::return_ptr)
    CALL return_ptr%Construct(name='soil_and_snow_properties', process_id=SSE_, owner_model_id=model_id)

  END FUNCTION Create_task_soil_and_snow_properties

  SUBROUTINE Register_sse_tasks(process, model_id)

    CLASS(t_jsb_process), INTENT(inout) :: process
    INTEGER,              INTENT(in)    :: model_id

    CALL process%Register_task(tsk_soil_and_snow_temp(model_id))
    CALL process%Register_task(tsk_soil_and_snow_properties(model_id))

  END SUBROUTINE Register_sse_tasks

  SUBROUTINE update_soil_and_snow_temperature(tile, options)

    USE mo_sse_process,            ONLY: calc_soil_temperature, calc_snow_temperature, calc_soiltemp_old
    USE mo_jsb_physical_constants, ONLY: alf, rhoh2o, dens_snow
    USE mo_sse_constants,          ONLY: snow_depth_min

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    !
    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(SEB_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: t_unfilt
    dsl4jsb_Real2D_onChunk :: snowmelt
    dsl4jsb_Real3D_onChunk :: t_soil_sl
    dsl4jsb_Real3D_onChunk :: t_soil_acoef
    dsl4jsb_Real3D_onChunk :: t_soil_bcoef
    dsl4jsb_Real3D_onChunk :: vol_heat_cap_sl
    dsl4jsb_Real3D_onChunk :: heat_cond_sl
    dsl4jsb_Real3D_onChunk :: t_snow
    dsl4jsb_Real3D_onChunk :: t_snow_acoef
    dsl4jsb_Real3D_onChunk :: t_snow_bcoef
    dsl4jsb_Real2D_onChunk :: vol_heat_cap_snow
    dsl4jsb_Real2D_onChunk :: heat_cond_snow
    dsl4jsb_Real2D_onChunk :: hcap_grnd
    dsl4jsb_Real2D_onChunk :: hcap_grnd_old
    dsl4jsb_Real2D_onChunk :: grnd_hflx
    dsl4jsb_Real2D_onChunk :: grnd_hflx_old
    dsl4jsb_Real3D_onChunk :: snow_depth_sl
    dsl4jsb_Real3D_onChunk :: soil_depth_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_sl
    dsl4jsb_Real3D_onChunk :: ice_soil_sl
    dsl4jsb_Real3D_onChunk :: wtr_freeze_sl
    dsl4jsb_Real3D_onChunk :: ice_melt_sl
    dsl4jsb_Real3D_onChunk :: wtr_soil_pot_scool_sl
    dsl4jsb_Real3D_onChunk :: vol_porosity_sl
    dsl4jsb_Real3D_onChunk :: matrix_pot_sl
    dsl4jsb_Real3D_onChunk :: bclapp_sl
    dsl4jsb_Real2D_onChunk :: thaw_depth
    dsl4jsb_Real2D_onChunk :: heat_cond
    dsl4jsb_Real2D_onChunk :: vol_heat_cap
    dsl4jsb_Real2D_onChunk :: weq_snow_soil
    dsl4jsb_Real2D_onChunk :: snow_soil_dens
    dsl4jsb_Real2D_onChunk :: fract_snow_soil

    ! Locally allocated vectors
    !
    REAL(wp), DIMENSION(options%nc) :: &
      & snow_depth,                    &
      & snow_fract,                    &
      & t_unfilt_corrected,            &
      & t_soil_top,                    &
      & hcap_grnd_snow,                &
      & hcap_grnd_soil,                &
      & grnd_hflx_snow,                &
      & grnd_hflx_soil
    LOGICAL, DIMENSION(options%nc)  :: &
      & ldglac
    REAL(wp), ALLOCATABLE :: ws_max_sl(:,:)

    INTEGER, DIMENSION(options%nc) ::  &
      & itop, itop_old

    LOGICAL  :: lstart
    INTEGER  :: iblk, ics, ice, nc, ic
    REAL(wp) :: dtime

    TYPE(t_jsb_vgrid), POINTER :: soil_e, snow_e
    INTEGER                    :: nsoil, is, nsnow, isnow

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_and_snow_temperature'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc
    dtime = options%dtime

    IF (.NOT. tile%Is_process_calculated(SSE_)) RETURN

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    model => Get_model(tile%owner_model_id)
    dsl4jsb_Get_config(SSE_)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      snow_e => Get_vgrid('snow_depth_energy')
      nsnow = snow_e%n_levels
    END IF

    ! Get reference to variables for current block
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(SEB_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var2D_onChunk(HYDRO_, snowmelt)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, weq_snow_soil)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)
    dsl4jsb_Get_var2D_onChunk(HYDRO_, fract_snow_soil)

    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_sl)         ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_acoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, t_soil_bcoef)      ! inout
    dsl4jsb_Get_var3D_onChunk(SSE_, vol_heat_cap_sl)   ! in
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_cond_sl)      ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap_snow) ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond_snow)    ! in
    dsl4jsb_Get_var2D_onChunk(SSE_, hcap_grnd)         ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, hcap_grnd_old)     ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, grnd_hflx)         ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, grnd_hflx_old)     ! out
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      dsl4jsb_Get_var3D_onChunk(SSE_, snow_depth_sl)   ! out
      dsl4jsb_Get_var3D_onChunk(SSE_, t_snow)          ! inout
      dsl4jsb_Get_var3D_onChunk(SSE_, t_snow_acoef)    ! inout
      dsl4jsb_Get_var3D_onChunk(SSE_, t_snow_bcoef)    ! inout
    END IF
    IF (.NOT. tile%is_glacier) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)  ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, matrix_pot_sl)    ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, bclapp_sl)        ! in
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)      ! inout
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_freeze_sl)    ! out
      dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_melt_sl)      ! out
      dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_pot_scool_sl) ! out
    END IF
    dsl4jsb_Get_var2D_onChunk(SSE_,   thaw_depth)      ! out

    dsl4jsb_Get_var2D_onChunk(SSE_,   heat_cond)
    dsl4jsb_Get_var2D_onChunk(SSE_,   vol_heat_cap)

    dsl4jsb_Get_var2D_onChunk(SEB_,   t_unfilt)        ! in

    lstart = is_time_experiment_start(options%current_datetime)

    !$ACC ENTER DATA CREATE(snow_depth, snow_fract, t_unfilt_corrected, t_soil_top) &
    !$ACC   CREATE(hcap_grnd_snow, hcap_grnd_soil, grnd_hflx_snow, grnd_hflx_soil, ldglac, itop, itop_old)

#ifdef _CRAYFTN
    ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
    !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
    DO ic = 1, nc
      ! Save old ground heat flux and heat capacity before updating
      ! The old flux and capacity (which have been used in the implicit calculation of the new surface temperature) are
      ! given back to the atmosphere so that the surface energy balance can be computed there
      grnd_hflx_old(ic) = grnd_hflx(ic)
      hcap_grnd_old(ic) = hcap_grnd(ic)

      ! Correction of surface temperature for snow and ice melt
      ! t_unfilt is later corrected in task "snowmelt_correction" of process SEB_, but we need it already here on the sub-tiles
      IF (lstart) THEN
        t_unfilt_corrected(ic) = t_unfilt(ic)
      ELSE
        t_unfilt_corrected(ic) = t_unfilt(ic) - snowmelt(ic) * dtime * alf / hcap_grnd(ic)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    IF (.NOT. dsl4jsb_Config(SSE_)%l_snow) THEN
      ldglac(:) = SPREAD(tile%is_glacier, DIM=1, ncopies=nc)
      !$ACC UPDATE DEVICE(ldglac) ASYNC(1)

      CALL calc_soiltemp_old(                          &
        & nc,                                          &
        & nsoil,                                       &
        & soil_e%dz(:),                                &
        & soil_e%levels(:),                            &
        & dtime,                                       &
        & lstart,                                      &
        & t_unfilt_corrected(:)                     ,  &
        & weq_snow_soil(:)                          ,  &
        & heat_cond(:)                              ,  &
        & vol_heat_cap(:)                           ,  &
        & t_soil_acoef                          (:,:), &
        & t_soil_bcoef                          (:,:), &
        & t_soil_sl                             (:,:), &
        & hcap_grnd                             (:),   &
        & grnd_hflx                             (:),   &
        & ldglac                                (:)    &
        !& SPREAD(tile%is_glacier, DIM=1, ncopies=nc)   &
        & )


      !$ACC WAIT(1)
      !$ACC EXIT DATA &
      !$ACC   DELETE(snow_depth, snow_fract, t_unfilt_corrected, t_soil_top) &
      !$ACC   DELETE(hcap_grnd_snow, hcap_grnd_soil, grnd_hflx_snow, grnd_hflx_soil, ldglac, itop, itop_old)
      RETURN

    END IF

    ! l_snow=.TRUE. below here

    IF (.NOT. tile%is_glacier) THEN

      ! Calculate snow depth in m (no longer in m water equivalent)
      IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
#ifdef _CRAYFTN
        ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
        !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
        DO ic = 1,nc
          snow_depth(ic) = weq_snow_soil(ic) * rhoh2o / snow_soil_dens(ic)
        END DO
        !$ACC END PARALLEL LOOP
      ELSE
#ifdef _CRAYFTN
        ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
        !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
        DO ic = 1,nc
          snow_depth(ic) = weq_snow_soil(ic) * rhoh2o / dens_snow
        END DO
        !$ACC END PARALLEL LOOP
      END IF

      ! Calculate top snow layer for previous time step
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
      DO ic = 1,nc
        itop_old(ic) = nsnow + 1 - COUNT(snow_depth_sl(ic,:) > 0._wp)
      END DO
      !$ACC END PARALLEL LOOP

      ! Calculate depths of snow layers
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1) COLLAPSE(2)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
#endif
      DO isnow=1,nsnow
        DO ic = 1,nc
          snow_depth_sl(ic,isnow) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
      DO ic = 1,nc
        snow_fract(ic) = fract_snow_soil(ic)
        IF (snow_depth(ic) > snow_depth_min .AND. snow_fract(ic) > EPSILON(1._wp)) THEN
          snow_depth(ic) = snow_depth(ic) / snow_fract(ic)
        ELSE
          snow_depth(ic) = 0._wp
          snow_fract(ic) = 0._wp
        END IF
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC WAIT(1)

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL PRESENT(sse__mem, seb__mem, hydro__mem) ASYNC(1)
#else
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
#endif
      !$ACC LOOP SEQ
      DO isnow=nsnow,1,-1
        !$ACC LOOP GANG VECTOR
        DO ic = 1,nc
          IF (snow_depth(ic) > 0._wp) THEN
            snow_depth_sl(ic,isnow) = MIN(snow_depth(ic), snow_e%dz(isnow))
            snow_depth(ic) = snow_depth(ic) - snow_depth_sl(ic,isnow)
          END IF
        END DO
        !$ACC END LOOP
      END DO
      !$ACC END LOOP
      !$ACC END PARALLEL

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
      DO ic = 1,nc
        snow_depth_sl(ic,nsnow) = snow_depth_sl(ic,nsnow) + snow_depth(ic)
      END DO
      !$ACC END PARALLEL LOOP

      ! Calculate top snow layer for current time step

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
      DO ic = 1,nc
        itop(ic) = nsnow + 1 - COUNT(snow_depth_sl(ic,:) > 0._wp)
      END DO
      !$ACC END PARALLEL LOOP

      CALL calc_snow_temperature( &
        ! in
        & nc,                     &
        & nsnow,                  &
        & soil_e%dz(1:2),         &
        & dtime,                  &
        & lstart,                 &
        & snow_depth_sl    (:,:), &
        & itop_old           (:), &
        & itop               (:), &
        & t_unfilt_corrected (:), &
        & vol_heat_cap_snow  (:), &
        & heat_cond_snow     (:), &
        & vol_heat_cap_sl  (:,1), &
        & heat_cond_sl     (:,1), &
        & t_soil_sl        (:,:), &
        & t_soil_acoef     (:,1), &
        & t_soil_bcoef     (:,1), &
        ! inout
        & t_snow           (:,:), &
        & t_snow_acoef     (:,:), &
        & t_snow_bcoef     (:,:), &
        ! out
        & t_soil_top         (:), &
        & hcap_grnd_snow     (:), &
        & grnd_hflx_snow     (:)  &
        & )

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
      DO ic = 1, nc
        IF (snow_depth_sl(ic, nsnow) > 0._wp) THEN  ! At least one snow layer above soil
          t_soil_top(ic) = snow_fract(ic) * t_soil_top(ic) + (1._wp - snow_fract(ic)) * t_unfilt_corrected(ic)
        ELSE                                ! No snow
          t_soil_top(ic) = t_unfilt_corrected(ic)
        END IF
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1) COLLAPSE(2)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
#endif
      DO isnow = 1,nsnow
        DO ic = 1, nc
            snow_depth_sl(ic,isnow) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
      DO ic = 1, nc
        snow_depth(ic) = weq_snow_soil(ic) * rhoh2o / dens_snow
        snow_fract(ic) = 0._wp
        t_soil_top(ic) = t_unfilt_corrected(ic)
        grnd_hflx_snow(ic) = 0._wp
        hcap_grnd_snow(ic) = 0._wp
      END DO
      !$ACC END PARALLEL LOOP

    END IF

    ! @todo: move water/ice calculations to HYDRO_ process
    IF (tile%is_glacier) THEN
      CALL calc_soil_temperature( &
        ! in
        & nc                                       , &
        & nsoil                                    , &
        & soil_e%dz(:)                             , &
        & soil_e%levels(:)                         , &
        & dtime                                    , &
        & lstart                                   , &
        & dsl4jsb_Config(SSE_)%l_freeze            , &
        & dsl4jsb_Config(SSE_)%l_supercool         , &
        & t_soil_top                          (:  ), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        ! inout
        & t_soil_sl                           (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd_soil                      (:  ), &
        & grnd_hflx_soil                      (:  ), &
        & thaw_depth                          (:  )  &
        & )
    ELSE
      ALLOCATE(ws_max_sl(nc,nsoil))
      !$ACC ENTER DATA CREATE(ws_max_sl)
#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1) COLLAPSE(2)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) COLLAPSE(2)
#endif
      DO is = 1,nsoil
        DO ic = 1, nc
          ws_max_sl(ic,is) = soil_depth_sl(ic,is) * vol_porosity_sl(ic,is)
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      CALL calc_soil_temperature( &
        ! in
        & nc                                       , &
        & nsoil                                    , &
        & soil_e%dz(:)                             , &
        & soil_e%levels(:)                         , &
        & dtime                                    , &
        & lstart                                   , &
        & dsl4jsb_Config(SSE_)%l_freeze            , &
        & dsl4jsb_Config(SSE_)%l_supercool         , &
        & t_soil_top                          (:  ), &
        & vol_heat_cap_sl                     (:,:), &
        & heat_cond_sl                        (:,:), &
        ! inout
        & t_soil_sl                           (:,:), &
        & t_soil_acoef                        (:,:), &
        & t_soil_bcoef                        (:,:), &
        ! out
        & hcap_grnd_soil                      (:  ), &
        & grnd_hflx_soil                      (:  ), &
        & thaw_depth                          (:  ), &
        ! optional, in
        & ws_max_sl                           (:,:), &
        & matrix_pot_sl                       (:,:), &
        & bclapp_sl                           (:,:), &
        ! optional, inout
        & wtr_soil_sl                         (:,:), &
        & ice_soil_sl                         (:,:), &
        ! optional, out
        & wtr_freeze_sl                       (:,:), &
        & ice_melt_sl                         (:,:), &
        & wtr_soil_pot_scool_sl               (:,:)  &
        & )
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(ws_max_sl)
      DEALLOCATE(ws_max_sl)
    END IF

#ifdef _CRAYFTN
      ! ACCWA (Cray Fortran 15.0.1) : present fails for certain types of pointers
      !$ACC PARALLEL LOOP PRESENT(sse__mem, seb__mem, hydro__mem) GANG VECTOR ASYNC(1)
#else
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
#endif
    DO ic = 1,nc
      IF (snow_depth_sl(ic,nsnow) > 0._wp) THEN
        grnd_hflx(ic) = snow_fract(ic) * grnd_hflx_snow(ic) + (1._wp - snow_fract(ic)) * grnd_hflx_soil(ic)
        hcap_grnd(ic) = snow_fract(ic) * hcap_grnd_snow(ic) + (1._wp - snow_fract(ic)) * hcap_grnd_soil(ic)
      ELSE
        grnd_hflx(ic) = grnd_hflx_soil(ic)
        hcap_grnd(ic) = hcap_grnd_soil(ic)
      END IF
    END DO
    !$ACC END PARALLEL LOOP

    !$ACC WAIT(1)
    !$ACC EXIT DATA DELETE(snow_depth, snow_fract, t_unfilt_corrected, t_soil_top, hcap_grnd_snow) &
    !$ACC   DELETE(hcap_grnd_soil, grnd_hflx_snow, grnd_hflx_soil, ldglac, itop, itop_old)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_and_snow_temperature

  SUBROUTINE aggregate_soil_and_snow_temperature(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    TYPE(t_jsb_model), POINTER :: model
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_memory(HYDRO_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_temperature'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_acoef,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_bcoef,     weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   t_soil_sl,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   grnd_hflx,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   grnd_hflx_old,    weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   hcap_grnd,        weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   hcap_grnd_old,    weighted_by_fract)
    IF (dsl4jsb_Config(SSE_)%l_snow) THEN
      dsl4jsb_Aggregate_onChunk(SSE_,   snow_depth_sl,    weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_,   t_snow_acoef,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_,   t_snow_bcoef,     weighted_by_fract)
      dsl4jsb_Aggregate_onChunk(SSE_,   t_snow,           weighted_by_fract)
    END IF
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_freeze_sl,         weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, ice_melt_sl,           weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(HYDRO_, wtr_soil_pot_scool_sl, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   thaw_depth,            weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   vol_heat_cap,          weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_,   heat_cond,             weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_temperature

  SUBROUTINE update_soil_and_snow_properties(tile, options)

    USE mo_jsb_physical_constants, ONLY: ci, dens_snow_min ! , cs

    USE mo_sse_constants,          ONLY: vol_hcap_ice, vol_hcap_snow, vol_hcap_org_top, vol_hcap_org_below, &
                                       & hcond_ice, hcond_snow, hcond_snow_min, &
                                       & hcond_org_top, hcond_dry_org_top, hcond_org_below, hcond_dry_org_below !, vol_hcap_snow_min
    USE mo_hydro_constants,        ONLY: vol_porosity_org_top, vol_porosity_org_below
    USE mo_sse_process,            ONLY: calc_vol_heat_capacity, calc_thermal_conductivity

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    ! Local variables
    dsl4jsb_Def_config(SSE_)
    dsl4jsb_Def_memory(SSE_)
    dsl4jsb_Def_config(HYDRO_)
    dsl4jsb_Def_memory(HYDRO_)

    ! Pointers to variables in memory
    dsl4jsb_Real2D_onChunk :: &
      & snow_soil_dens,       &
      & vol_heat_cap,         &
      & heat_cond,            &
      & vol_heat_cap_snow,    &
      & heat_cond_snow
    dsl4jsb_Real3D_onChunk :: &
      & soil_depth_sl,        &
      & fract_org_sl,         &
      & vol_heat_cap_sl,      &
      & heat_cond_sl,         &
      & wtr_soil_sl,          &
      & ice_soil_sl,          &
      & vol_porosity_sl,      &
      & vol_field_cap_sl

    ! Locally allocated vectors
    !
    dsl4jsb_Real3D_onChunk :: &
      & fract_water,          &
      & fract_ice

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_vgrid), POINTER :: soil_e

    REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: &
      & vol_heat_cap_tmp,                    &
      & vol_heat_cap_org_tmp,                &
      & vol_porosity_org_tmp,                &
      & heat_cond_tmp,                       &
      & hcond_org_tmp,                       &
      & hcond_dry_org_tmp

    INTEGER  :: iblk, ics, ice, nc, nsoil, ic, is
    REAL(wp), POINTER :: soil_e_dz(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_soil_and_snow_properties'

    iblk  = options%iblk
    ics   = options%ics
    ice   = options%ice
    nc    = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    model => Get_model(tile%owner_model_id)

    soil_e => Get_vgrid('soil_depth_energy')
    nsoil = soil_e%n_levels
    soil_e_dz => soil_e%dz

    dsl4jsb_Get_config(SSE_)
    dsl4jsb_Get_config(HYDRO_)

    ! Get reference to variables for domain
    !
    dsl4jsb_Get_memory(SSE_)
    dsl4jsb_Get_memory(HYDRO_)

    dsl4jsb_Get_var3D_onChunk(SSE_, vol_heat_cap_sl)   ! out
    dsl4jsb_Get_var3D_onChunk(SSE_, heat_cond_sl)      ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap_snow) ! out
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond_snow) ! out


    ! On glacier tile, set variables and return
    IF (tile%is_glacier) THEN
      IF (model%config%l_compat401) THEN
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO is = 1, nsoil
          DO ic = 1, nc
            vol_heat_cap_sl(ic,is) = 2.09e+6_wp
            heat_cond_sl   (ic,is) = 2.508_wp
          END DO
        END DO
        !$ACC END LOOP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          vol_heat_cap_snow(ic) = 634500._wp
          heat_cond_snow   (ic) = 0.31_wp
        END DO
        !$ACC END LOOP
        !$ACC END PARALLEL
      ELSE
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1)
        !$ACC LOOP GANG(STATIC: 1) VECTOR COLLAPSE(2)
        DO is = 1, nsoil
          DO ic = 1, nc
            vol_heat_cap_sl(ic,is) = vol_hcap_ice
            heat_cond_sl   (ic,is) = hcond_ice
          END DO
        END DO
        !$ACC END LOOP
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO ic = 1, nc
          vol_heat_cap_snow(ic) = vol_hcap_snow
          heat_cond_snow   (ic) = hcond_snow
        END DO
        !$ACC END LOOP
        !$ACC END PARALLEL
      END IF
      IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')
      RETURN
    END IF

    dsl4jsb_Get_var3D_onChunk(HYDRO_, soil_depth_sl)
    dsl4jsb_Get_var2D_onChunk(SSE_, vol_heat_cap)
    dsl4jsb_Get_var2D_onChunk(SSE_, heat_cond)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, wtr_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, ice_soil_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_porosity_sl)
    dsl4jsb_Get_var3D_onChunk(HYDRO_, vol_field_cap_sl)

    IF (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil) THEN
      dsl4jsb_Get_var3D_onChunk(HYDRO_, fract_org_sl)
    ELSE
      ALLOCATE(fract_org_sl(nc, nsoil))
      !$ACC ENTER DATA CREATE(fract_org_sl)
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          fract_org_sl(ic,is) = 0._wp
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (dsl4jsb_Config(SSE_)%l_heat_cap_dyn) THEN
      ALLOCATE(fract_water(nc, nsoil), fract_ice(nc, nsoil))
      !$ACC ENTER DATA CREATE(fract_water, fract_ice)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          IF (soil_depth_sl(ic,is) > 0._wp) THEN
            fract_water(ic,is) = wtr_soil_sl(ic,is) / soil_depth_sl(ic,is)
            fract_ice  (ic,is) = ice_soil_sl(ic,is) / soil_depth_sl(ic,is)
          ELSE
            fract_water(ic,is) = 0._wp
            fract_ice  (ic,is) = 0._wp
          END IF
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      ALLOCATE(vol_heat_cap_tmp(nc,nsoil), vol_heat_cap_org_tmp(nc,nsoil))
      !$ACC ENTER DATA CREATE(vol_heat_cap_tmp, vol_heat_cap_org_tmp)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        vol_heat_cap_tmp    (ic,1) = vol_heat_cap(ic)
        vol_heat_cap_org_tmp(ic,1) = vol_hcap_org_top
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 2, nsoil
        DO ic = 1, nc
          vol_heat_cap_tmp    (ic,is) = vol_heat_cap(ic)
          vol_heat_cap_org_tmp(ic,is) = vol_hcap_org_below
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          vol_heat_cap_sl(ic,is) = calc_vol_heat_capacity(          &
            &                         vol_heat_cap_tmp(ic,is)     , &
            &                         vol_heat_cap_org_tmp(ic,is) , &
            &                         vol_porosity_sl(ic,is)      , &
            &                         fract_org_sl(ic,is)         , &
            &                         fract_water(ic,is)          , &
            &                         fract_ice(ic,is)            , &
            &                         soil_depth_sl(ic,is)        , &
            &                         soil_e_dz(is)                 &
            & )
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(fract_water, fract_ice, vol_heat_cap_tmp, vol_heat_cap_org_tmp)
      DEALLOCATE(fract_water, fract_ice, vol_heat_cap_tmp, vol_heat_cap_org_tmp)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          vol_heat_cap_sl(ic,is) = vol_heat_cap(ic)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    ENDIF

    IF (dsl4jsb_Config(SSE_)%l_heat_cond_dyn) THEN

      ALLOCATE(heat_cond_tmp(nc,nsoil), vol_porosity_org_tmp(nc,nsoil), &
             & hcond_org_tmp(nc,nsoil), hcond_dry_org_tmp(nc,nsoil))
      !$ACC ENTER DATA CREATE(heat_cond_tmp, vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp)

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        heat_cond_tmp(ic,1)         = heat_cond(ic)
        vol_porosity_org_tmp(ic,1)  = vol_porosity_org_top
        hcond_org_tmp(ic,1)         = hcond_org_top
        hcond_dry_org_tmp(ic,1)     = hcond_dry_org_top
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 2, nsoil
        DO ic = 1, nc
          heat_cond_tmp(ic,is)        = heat_cond(ic)
          vol_porosity_org_tmp(ic,is) = vol_porosity_org_below
          hcond_org_tmp(ic,is)        = hcond_org_below
          hcond_dry_org_tmp(ic,is)    = hcond_dry_org_below
        END DO
      END DO
      !$ACC END PARALLEL LOOP

      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          heat_cond_sl (ic,is) = calc_thermal_conductivity( &
            &                                heat_cond_tmp(ic,is)          , &
            &                                wtr_soil_sl(ic,is)            , &
            &                                ice_soil_sl(ic,is)            , &
            &                                vol_porosity_sl(ic,is)        , &
            &                                vol_field_cap_sl(ic,is)       , &
            &                                fract_org_sl(ic,is)           , &
            &                                vol_porosity_org_tmp(ic,is)   , &
            &                                hcond_org_tmp(ic,is)          , &
            &                                hcond_dry_org_tmp(ic,is)      , &
            &                                soil_depth_sl(ic,is)          , &
            &                                soil_e_dz(is)                   &
            & )
        END DO
      END DO
      !$ACC END PARALLEL LOOP
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(heat_cond_tmp, vol_porosity_org_tmp, hcond_org_tmp, hcond_dry_org_tmp)
      DEALLOCATE(heat_cond_tmp, vol_porosity_org_tmp,hcond_org_tmp, hcond_dry_org_tmp)
    ELSE
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR COLLAPSE(2) ASYNC(1)
      DO is = 1, nsoil
        DO ic = 1, nc
          heat_cond_sl(ic,is) = heat_cond(ic)
        END DO
      END DO
      !$ACC END PARALLEL LOOP
    END IF

    IF (.NOT. (dsl4jsb_Config(HYDRO_)%l_organic .AND. tile%contains_soil)) THEN
      !$ACC WAIT(1)
      !$ACC EXIT DATA DELETE(fract_org_sl)
      DEALLOCATE(fract_org_sl)
    END IF

    IF (model%config%l_compat401) THEN
      !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
      DO ic = 1, nc
        vol_heat_cap_snow(ic) = 634500._wp
        heat_cond_snow(ic) = 0.31_wp
      END DO
      !$ACC END PARALLEL LOOP
    ELSE
      IF (dsl4jsb_Config(SSE_)%l_snow) THEN
        IF (dsl4jsb_Config(SSE_)%l_dynsnow) THEN
          dsl4jsb_Get_var2D_onChunk(HYDRO_, snow_soil_dens)
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO ic = 1, nc
            IF (snow_soil_dens(ic) > dens_snow_min) THEN         ! "old" or mixture of old and fresh snow
              vol_heat_cap_snow(ic) = ci * snow_soil_dens(ic)    ! Use specific heat of ice (Verseghy 1991)
              ! heat_cond_snow   (:) = 2.576e-6_wp * snow_soil_dens(:) ** 2._wp + 0.074_wp                      ! Verseghy (1991)
              ! heat_cond_snow   (:) = 2.9e-6_wp*snow_soil_dens(:)**2._wp                                       ! Goodrich (1982)
              ! Calonne et al. (2011):
              heat_cond_snow   (ic) = 2.5e-6_wp*snow_soil_dens(ic)**2._wp - 0.123e-3_wp*snow_soil_dens(ic) + 0.024_wp
            ELSE                                         ! fresh snow
              ! @todo: these minimum values from mo_sse_constants.f90 don't match with the above equations and minimum snow density
              ! vol_heat_cap_snow(:) = vol_hcap_snow_min
              heat_cond_snow   (ic) = hcond_snow_min
              vol_heat_cap_snow(ic) = ci * dens_snow_min
            END IF
          END DO
          !$ACC END PARALLEL LOOP
        ELSE
          !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
          DO ic = 1, nc
            vol_heat_cap_snow(ic) = vol_hcap_snow  ! Use specific heat of snow (cs * dens_snow)
            heat_cond_snow   (ic) = hcond_snow
          END DO
          !$ACC END PARALLEL LOOP
        END IF
      ELSE
        !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
        DO ic = 1, nc
          vol_heat_cap_snow(ic) = vol_hcap_snow    ! Use specific heat of snow (cs * dens_snow)
          heat_cond_snow   (ic) = hcond_snow
        END DO
        !$ACC END PARALLEL LOOP
      END IF
    END IF

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE update_soil_and_snow_properties
  !
  ! -------------------------------------------------------------------------------------------------------
  !>
  !! Implementation of "aggregate" for task "soil_properties"
  !!
  !! @param[in,out] tile    Tile for which routine is executed.
  !! @param[in]     options Additional run-time parameters.
  !!
  SUBROUTINE aggregate_soil_and_snow_properties(tile, options)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options

    dsl4jsb_Def_memory(SSE_)

    CLASS(t_jsb_aggregator), POINTER :: weighted_by_fract

    CHARACTER(len=*), PARAMETER :: routine = modname//':aggregate_soil_and_snow_properties'

    INTEGER :: iblk, ics, ice

    iblk = options%iblk
    ics  = options%ics
    ice  = options%ice

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(SSE_)

    weighted_by_fract => tile%Get_aggregator("weighted_by_fract")

    dsl4jsb_Aggregate_onChunk(SSE_, vol_heat_cap_sl,   weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, heat_cond_sl,      weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, vol_heat_cap_snow, weighted_by_fract)
    dsl4jsb_Aggregate_onChunk(SSE_, heat_cond_snow,    weighted_by_fract)

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE aggregate_soil_and_snow_properties

#endif
END MODULE mo_sse_interface
