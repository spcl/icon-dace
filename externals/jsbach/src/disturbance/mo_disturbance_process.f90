!> Contains the routines for the disturb processes
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
MODULE mo_disturb_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: burned_fract_jsbach, broken_woody_fract_jsbach, get_relative_humidity_air

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_process'

CONTAINS

  SUBROUTINE burned_fract_jsbach( &
    & fire_rel_hum_threshold,                    & ! in
    & fire_litter_threshold,                     & ! in
    & fire_minimum,                              & ! in
    & fire_tau,                                  & ! in
    & q_rel_air_climbuf,                         & ! in
    & fuel,                                      & ! in
    & burned_fract                               & ! inout
    & )

    ! Input Arguments
    REAL(wp), INTENT(in) ::     &
      & fire_rel_hum_threshold, &  ! Maximum relative humidity for fire
      & fire_litter_threshold,  &  ! Minimum amount of litter necessary for fire [mol(C)/m^2(grid box)]
      & fire_minimum,           &  ! Minimum fraction woody PFTs burning every year
      & fire_tau,               &  ! Return period of fire at 0% relative humidity [year]
      & q_rel_air_climbuf(:),   &  ! Relative humidity (smoothed) [%]
      & fuel(:)                    ! Amount of fuel [mol(C)/m^2(grid box)]

    ! Output Arguments
    REAL(wp), INTENT(inout) :: burned_fract(:)

    ! Locals
    INTEGER  :: ic, nc

    REAL(wp), PARAMETER :: delta_time_yr = 1._wp / 365._wp

    nc = SIZE(fuel)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      ! Minimum fraction burnt every day (needed with natural LCC)
      burned_fract(ic) = fire_minimum * delta_time_yr

      ! If conditions for fire are given ...
      IF (q_rel_air_climbuf(ic) < fire_rel_hum_threshold .AND. fuel(ic) > fire_litter_threshold) THEN
        burned_fract(ic) = burned_fract(ic)                                                                        &
          &               + (fire_rel_hum_threshold - q_rel_air_climbuf(ic)) / (fire_tau * fire_rel_hum_threshold) &
          &                 * delta_time_yr
      END IF

    END DO
    !$END PARALLEL LOOP

  END SUBROUTINE burned_fract_jsbach

  ! Calculate windbreak for woody types
  SUBROUTINE broken_woody_fract_jsbach( &
    & wind_threshold,                                  & ! in
    & wind_damage_scale,                               & ! in
    & cover_fract_pot,                                 & ! in
    & prev_day_max_wind_10m,                           & ! in
    & max_wind_10m,                                    & ! in
    & damaged_fract                                    & ! inout
    & )

    ! Arguments
    REAL(wp), INTENT(in) ::       &
      & wind_threshold,           &
      & wind_damage_scale,        &
      & cover_fract_pot(:),       &
      & prev_day_max_wind_10m(:), &
      & max_wind_10m(:)

    REAL(wp), INTENT(inout) :: &
      & damaged_fract(:)

    ! Locals
    INTEGER  :: ic, nc

    REAL(wp), PARAMETER :: delta_time_yr = 1._wp / 365._wp

    nc = SIZE(max_wind_10m)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      IF (cover_fract_pot(ic) > EPSILON(1._wp) .AND. prev_day_max_wind_10m(ic) > max_wind_10m(ic) * wind_threshold) THEN
        damaged_fract(ic) = MIN( 1._wp - EPSILON(1._wp) / cover_fract_pot(ic),                                          &
        &                        delta_time_yr * wind_damage_scale * prev_day_max_wind_10m(ic) ** 3._wp / max_wind_10m(ic) )
      ENDIF
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE broken_woody_fract_jsbach

  FUNCTION get_relative_humidity_air(air_moisture, air_temperature, air_press) RESULT(relative_humidity_air)

    USE mo_jsb_physical_constants, ONLY: rd  ! gas constant for dry air [J/(K*kg)]
    USE mo_jsb_physical_constants, ONLY: rv  ! gas constant for water vapor [J/(K*kg)]
    USE mo_phy_schemes,            ONLY: qsat_water

    REAL(wp), INTENT(in) :: &
      air_moisture(:),      & ! Specific humidity of air [kg kg-1]
      air_temperature(:),   & ! Temperature of air [K]
      air_press(:)            ! Pressure of air [Pa]

    REAL(wp) :: air_qsat      ! Saturation specific humidity
    REAL(wp) :: relative_humidity_air(SIZE(air_moisture))

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_relative_humidity_air'
    INTEGER :: ic, nc

    nc = SIZE(air_moisture)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1) &
    !$ACC   PRIVATE(air_qsat)
    DO ic = 1, nc
      air_qsat = qsat_water(air_temperature(ic), air_press(ic), use_convect_tables=.TRUE.)
      relative_humidity_air(ic) = 100._wp * air_moisture(ic)                                       &
                                 * (  (1._wp - rd / rv) * air_qsat + (rd / rv) )                   &
                                 / ( ((1._wp - rd / rv) * air_moisture(ic) + (rd / rv)) * air_qsat )
      relative_humidity_air(ic) = MIN(100._wp, relative_humidity_air(ic))
    END DO
    !$ACC END PARALLEL LOOP

  END FUNCTION get_relative_humidity_air

#endif
END MODULE mo_disturb_process
