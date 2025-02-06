!> lper routines for assimilation
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
!>#### various helper routines for the assimilation process
!>
MODULE mo_q_assimi_util
#ifndef __NO_QUINCY__

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reset_q_assimi_fluxes, calculate_time_average_q_assimilation

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_assimi_util'


CONTAINS


  ! ====================================================================================================== !
  !
  !> init/reset assimilation fluxes (with/to zero) - called prior to other assimi tasks
  !
  SUBROUTINE reset_q_assimi_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_math_constants,    ONLY: zero

    dsl4jsb_Use_processes Q_ASSIMI_
    dsl4jsb_Use_memory(Q_ASSIMI_)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_q_assimi_fluxes'

    dsl4jsb_Def_memory(Q_ASSIMI_)
    ! -------------------------------------------------------------------------------------------------- !
    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN

    model  => Get_model(tile%owner_model_id)

    IF (model%lctlib(tile%lcts(1)%lib_id)%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(Q_ASSIMI_)

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! set zero all flux variables of the assimilation process
    dsl4jsb_var2D_onChunk(Q_ASSIMI_, gross_assimilation)           = zero
    dsl4jsb_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C13)       = zero
    dsl4jsb_var2D_onChunk(Q_ASSIMI_, gross_assimilation_C14)       = zero
    dsl4jsb_var2D_onChunk(Q_ASSIMI_, net_assimilation)             = zero
    dsl4jsb_var2D_onChunk(Q_ASSIMI_, net_assimilation_boc)         = zero
    dsl4jsb_var2D_onChunk(Q_ASSIMI_, maint_respiration_leaf)       = zero

    dsl4jsb_var3D_onChunk(Q_ASSIMI_, gross_assimilation_cl)      = zero
    dsl4jsb_var3D_onChunk(Q_ASSIMI_, net_assimilation_cl)        = zero
    dsl4jsb_var3D_onChunk(Q_ASSIMI_, maint_respiration_leaf_cl)  = zero

    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE reset_q_assimi_fluxes


  ! ======================================================================================================= !
  !> calculate time moving averages and daytime averages for Q_ASSIMI_
  !>
  SUBROUTINE calculate_time_average_q_assimilation(tile, options)
    USE mo_kind,                    ONLY: wp
    USE mo_jsb_math_constants,      ONLY: eps8, eps1
    USE mo_jsb_physical_constants,  ONLY: Tzero
    USE mo_jsb_impl_constants,      ONLY: test_false_true
    USE mo_jsb_control,             ONLY: debug_on
    USE mo_jsb_class,               ONLY: Get_model
    USE mo_jsb_model_class,         ONLY: t_jsb_model
    USE mo_jsb_tile_class,          ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,          ONLY: t_jsb_task_options
    USE mo_lnd_time_averages        ! e.g. calc_time_mavg, mavg_period_tphen, mavg_period_weekly
    dsl4jsb_Use_processes A2L_, Q_ASSIMI_, Q_PHENO_
    dsl4jsb_Use_memory(A2L_)
    dsl4jsb_Use_memory(Q_ASSIMI_)
    dsl4jsb_Use_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    TYPE(t_jsb_task_options),   INTENT(in)    :: options
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER              :: model
    LOGICAL,  DIMENSION(options%nc)         :: l_growing_season     !< growing_season LOGICAL
    REAl(wp), DIMENSION(options%nc)         :: hlp                  !< helper
    REAL(wp)                                :: dtime                !< timestep
    INTEGER                                 :: ic                   !< looping
    INTEGER                                 :: iblk, ics, ice, nc   !< dimensions
    CHARACTER(len=*), PARAMETER             :: routine = modname//':calculate_time_average_q_assimilation'
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Def_memory(A2L_)
    dsl4jsb_Def_memory(Q_ASSIMI_)
    dsl4jsb_Def_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Real2D_onChunk      :: t_air
    dsl4jsb_Real2D_onChunk      :: swpar_srf_down
    dsl4jsb_Real2D_onChunk      :: daytime_counter
    dsl4jsb_Real2D_onChunk      :: local_time_day_seconds
    ! Q_ASSIMI_
    dsl4jsb_Real2D_onChunk      :: beta_air
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs
    dsl4jsb_Real2D_onChunk      :: beta_soa
    dsl4jsb_Real2D_onChunk      :: soa_tsoa_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soa_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tphen_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_tfrac_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_tcnl_mavg
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime
    dsl4jsb_Real2D_onChunk      :: beta_air_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_ps_daytime_dacc
    dsl4jsb_Real2D_onChunk      :: beta_soil_gs_daytime_dacc
    ! Q_PHENO_
    dsl4jsb_Real2D_onChunk      :: growing_season
    ! ----------------------------------------------------------------------------------------------------- !
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc
    dtime   = options%dtime
    ! ----------------------------------------------------------------------------------------------------- !
    IF (.NOT. tile%Is_process_calculated(Q_ASSIMI_)) RETURN
    IF (tile%lcts(1)%lib_id == 0) RETURN !< only if the present tile is a pft
    ! ----------------------------------------------------------------------------------------------------- !
    model  => Get_model(tile%owner_model_id)
    ! ----------------------------------------------------------------------------------------------------- !
    IF (model%lctlib(tile%lcts(1)%lib_id)%BareSoilFlag) RETURN !< do not run this routine at tiles like "bare soil" and "urban area"
    IF (debug_on() .AND. iblk == 1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')
    ! ----------------------------------------------------------------------------------------------------- !
    dsl4jsb_Get_memory(A2L_)
    dsl4jsb_Get_memory(Q_ASSIMI_)
    dsl4jsb_Get_memory(Q_PHENO_)
    ! ----------------------------------------------------------------------------------------------------- !
    ! A2L_
    dsl4jsb_Get_var2D_onChunk(A2L_, t_air)                            ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, swpar_srf_down)                   ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, daytime_counter)                  ! in
    dsl4jsb_Get_var2D_onChunk(A2L_, local_time_day_seconds)           ! in
    ! Q_ASSIMI_
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air)                    ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps)                ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs)                ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soa)                    ! in
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, soa_tsoa_mavg)               ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soa_tphen_mavg)         ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tphen_mavg)     ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps_tfrac_mavg)     ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tfrac_mavg)     ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps_tcnl_mavg)      ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_tcnl_mavg)      ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air_tfrac_mavg)         ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air_tcnl_mavg)          ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air_daytime)            ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps_daytime)        ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_daytime)        ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_air_daytime_dacc)       ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_ps_daytime_dacc)   ! out
    dsl4jsb_Get_var2D_onChunk(Q_ASSIMI_, beta_soil_gs_daytime_dacc)   ! out
    ! Q_PHENO_
    dsl4jsb_Get_var2D_onChunk(Q_PHENO_, growing_season)               ! in


    !>0.9 transform REAL growing_season(:) into LOGICAL l_growing_season(:)/_cl(:,:)
    !>
    DO ic = 1, nc
      IF (growing_season(ic) > test_false_true) THEN
        l_growing_season(ic) = .TRUE.
      ELSE
        l_growing_season(ic) = .FALSE.
      END IF
    END DO

    ! ----------------------------------------------------------------------------------------------------- !
    !>1.0 daytime averages
    !>
    !>

    !>  1.1 accumulate values - if light is available (i.e. daytime)
    !>
    DO ic = 1, nc
      IF (swpar_srf_down(ic) > eps8) THEN
        ! add values
        beta_air_daytime_dacc(ic)       = beta_air_daytime_dacc(ic)     + beta_air(ic)
        beta_soil_ps_daytime_dacc(ic)   = beta_soil_ps_daytime_dacc(ic) + beta_soil_ps(ic)
        beta_soil_gs_daytime_dacc(ic)   = beta_soil_gs_daytime_dacc(ic) + beta_soil_gs(ic)
      END IF
    END DO

    !>  1.2 calc daytime averages - at the timestep after local midnight (1st timestep of the new day)
    !>
    !>    calculate the average of the previous day and pass this value to the according variable
    DO ic = 1, nc
      IF (ABS(local_time_day_seconds(ic) - dtime) < eps1) THEN
        ! check if at least one value had been assigned at the current day, avoiding division by zero
        IF (daytime_counter(ic) > eps8) THEN
          beta_air_daytime(ic)     = beta_air_daytime_dacc(ic)     / daytime_counter(ic)
          beta_soil_ps_daytime(ic) = beta_soil_ps_daytime_dacc(ic) / daytime_counter(ic)
          beta_soil_gs_daytime(ic) = beta_soil_gs_daytime_dacc(ic) / daytime_counter(ic)
        END IF
        ! zero the accumulation variables after daily average has been calculated
        beta_air_daytime_dacc(ic)     = 0.0_wp
        beta_soil_ps_daytime_dacc(ic) = 0.0_wp
        beta_soil_gs_daytime_dacc(ic) = 0.0_wp
      END IF
    END DO

    ! ----------------------------------------------------------------------------------------------------- !
    !>2.0 moving averages
    !>
    !>

    ! docu:
    ! calc_time_mavg(dtime, current average, new value, length of avg_period,  !
    !                do_calc=LOGICAL, avg_period_unit='day')            ! OPTIONAL
    !                RETURN(new current average)
    ! the unit of the averaging period is 'day' by default, but can also be 'week' or 'year'

    !>  2.1 tphen (averages at the time-scale of phenology)
    !>
    beta_soil_gs_tphen_mavg(:) = calc_time_mavg(dtime, beta_soil_gs_tphen_mavg(:), beta_soil_gs(:), mavg_period_tphen)
    beta_soa_tphen_mavg(:)     = calc_time_mavg(dtime, beta_soa_tphen_mavg(:), beta_soa(:), mavg_period_tphen)

    !>  2.2 soa time-scale
    !>
    hlp(:)                     = t_air(:) - Tzero     ! convert from Kelvin to deg Celsius
    soa_tsoa_mavg(:)           = calc_time_mavg(dtime, soa_tsoa_mavg(:), hlp(:), mavg_period_tsoa)

    !>  2.3 tfrac (averages at the time-scale of within-leaf N allocation fractions)
    !>
    beta_soil_ps_tfrac_mavg(:) = calc_time_mavg(dtime, beta_soil_ps_tfrac_mavg(:), beta_soil_ps_daytime(:), &
                                                mavg_period_tfrac)
    beta_soil_gs_tfrac_mavg(:) = calc_time_mavg(dtime, beta_soil_gs_tfrac_mavg(:), beta_soil_gs_daytime(:), &
                                                mavg_period_tfrac)
    beta_air_tfrac_mavg(:)     = calc_time_mavg(dtime, beta_air_tfrac_mavg(:), beta_air_daytime(:), &
                                                mavg_period_tfrac)

    !>  2.4 tcnl (averages at the time-scale of leaf N allocation fractions)
    !>
    beta_soil_ps_tcnl_mavg(:)  = calc_time_mavg(dtime, beta_soil_ps_tcnl_mavg(:), beta_soil_ps(:), &
                                                mavg_period_tcnl, do_calc=l_growing_season(:))
    beta_soil_gs_tcnl_mavg(:)  = calc_time_mavg(dtime, beta_soil_gs_tcnl_mavg(:), beta_soil_gs(:), &
                                                mavg_period_tcnl, do_calc=l_growing_season(:))
    beta_air_tcnl_mavg(:)      = calc_time_mavg(dtime, beta_air_tcnl_mavg(:), beta_air(:), &
                                                mavg_period_tcnl, do_calc=l_growing_season(:))

  END SUBROUTINE calculate_time_average_q_assimilation

#endif
END MODULE mo_q_assimi_util
