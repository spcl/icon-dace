!> quincy atm2land routines: only required within ICON-Land and not in QS (thus _iq_)
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
!>### Contains routines for the a2l (atm2land) proc
!>
MODULE mo_iq_atm2land_process
#ifndef __NO_QUINCY__

  USE mo_kind,                  ONLY: wp
  USE mo_jsb_math_constants,    ONLY: eps1, eps8

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_local_time_and_daytime_counter, update_slow_sb_pool_accelerator_bookkeeping

  CHARACTER(len=*), PARAMETER :: modname = 'mo_iq_atm2land_process'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> updates local time and updates or resets quincys daytime counter depending on the local time
  !>
  PURE ELEMENTAL SUBROUTINE update_local_time_and_daytime_counter( &
      &                       global_seconds_day, dtime, lon, swpar_srf_down, daytime_counter, local_time_day_seconds)

    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER ,  INTENT(in)        :: global_seconds_day     !< seconds passed on this day (global time)
    REAL(wp),  INTENT(in)        :: dtime                  !< timestep
    REAL(wp),  INTENT(in)        :: lon                    !< longitude
    REAL(wp),  INTENT(in)        :: swpar_srf_down         !< current shortwave surface downward flux
    REAL(wp),  INTENT(inout)     :: daytime_counter        !< nr of timesteps of current day with sunlight
    REAL(wp),  INTENT(inout)     :: local_time_day_seconds !< seconds passed on this day (local cell time)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER  :: time_shift_in_hours

    CHARACTER(len=*), PARAMETER :: routine = modname//':update_local_time_and_daytime_counter'
    ! ----------------------------------------------------------------------------------------------------- !

    ! reset the daytime_counter if last timestep was midnight
    IF (ABS(local_time_day_seconds - dtime) < eps1) THEN
      daytime_counter = 0.0_wp
    ENDIF

    IF (swpar_srf_down > eps8) THEN
      ! count number of daytime values
      daytime_counter = daytime_counter + 1.0_wp
    ENDIF

    ! calculate the local time shift in whole hours (NINT!) relative to the global time (15 degrees per hour)
    time_shift_in_hours = NINT(MODULO(lon,360._wp)/15._wp)
    ! convert the time shift from whole hours to whole seconds and calculate the seconds passed on this day in local time
    ! (MODULO 86400 is applied to detemine the seconds of the local day, as this can be the day before or after the global day)
    local_time_day_seconds  = REAL(MODULO(global_seconds_day + (time_shift_in_hours * 3600), 86400),wp)

  END SUBROUTINE update_local_time_and_daytime_counter

  ! ======================================================================================================= !
  !>
  !> updates counters and flag for the sb slow pool spin-up accelerator
  !>
  PURE ELEMENTAL SUBROUTINE update_slow_sb_pool_accelerator_bookkeeping( dtime, current_year, max_executions, frequency, start_year, &
      &                                                                execution_counter, execute_accelerator)
    ! ----------------------------------------------------------------------------------------------------- !
    REAL(wp),  INTENT(in)    :: dtime               !< timestep
    INTEGER ,  INTENT(in)    :: current_year        !< simulation year
    INTEGER ,  INTENT(in)    :: max_executions      !< maximum number of executions of the accelerator
    INTEGER ,  INTENT(in)    :: frequency           !< frequency (years) in which the accelerator is executed
    INTEGER ,  INTENT(in)    :: start_year          !< years after which the spin-up accelerator should first be executed
    REAL(wp),  INTENT(inout) :: execution_counter   !< the number of timesteps since the beginning of spinup [#]
    REAL(wp),  INTENT(inout) :: execute_accelerator !< represents logical: execute the spin-up accelerator (0 = FALSE, 1= True)
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':update_slow_sb_pool_accelerator_bookkeeping'
    ! ----------------------------------------------------------------------------------------------------- !
    !> stop increasing the time counter if the maximum execution time is already exceeded
    IF(execution_counter <= max_executions) THEN
      !> if accelerator was executed, set the logical value 'execute' as True (1.0_wp); the 'execution counter' +1
      IF ((current_year >= start_year) .AND. (MOD(current_year, frequency) < eps1 )) THEN
        execute_accelerator = 1.0_wp
        execution_counter = execution_counter + 1.0_wp
      ENDIF
    ENDIF
  END SUBROUTINE update_slow_sb_pool_accelerator_bookkeeping


#endif
END MODULE mo_iq_atm2land_process

