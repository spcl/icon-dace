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
  USE mo_exception, ONLY: finish, message_text

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: burned_fract_jsbach, broken_woody_fract_jsbach, get_relative_humidity_air

  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_process'

CONTAINS

  PURE ELEMENTAL SUBROUTINE burned_fract_jsbach( &
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
      & fire_rel_hum_threshold, &
      & fire_litter_threshold,  &
      & fire_minimum,           &
      & fire_tau,               &
      & q_rel_air_climbuf,      &
      & fuel

    ! Output Arguments
    REAL(wp), INTENT(inout) :: burned_fract

    ! Locals
    REAL(wp) :: delta_time_yr
                !fire_rel_hum_threshold, & ! maximal relative humidity for fire
                !fire_litter_threshold,  & ! minimal amount of litter [mol(C)/m^2(grid box)] for fire
                !fire_minimum_woody,     & ! minimal fraction of act_fpc of woody PFT to be burned each year
                !fire_minimum_grass,     & ! minimal fraction of act_fpc of grass PFT to be burned each year
                !fire_tau_woody,         & ! return period of fire for woody PFT [year] at 0% relative humidity
                !fire_tau_grass,         & ! return period of fire for grass PFT [year] at 0% relative humidity

    delta_time_yr = 1._wp / 365._wp

    burned_fract = fire_minimum * delta_time_yr

    IF (       q_rel_air_climbuf < fire_rel_hum_threshold   &
         & .AND. fuel > fire_litter_threshold               &
       ) THEN
      burned_fract = burned_fract                                                                        &
        &           + (fire_rel_hum_threshold - q_rel_air_climbuf) / (fire_tau * fire_rel_hum_threshold) &
        &             * delta_time_yr
    END IF

  END SUBROUTINE burned_fract_jsbach

  ! Calculate windbreak for woody types
  PURE ELEMENTAL SUBROUTINE broken_woody_fract_jsbach( &
    & wind_threshold,                                  & ! in
    & wind_damage_scale,                               & ! in
    & cover_fract_pot,                                 & ! in
    & prev_day_max_wind_10m,                           & ! in
    & max_wind_10m,                                    & ! in
    & damaged_fract                                    & ! inout
    & )

    ! Arguments
    REAL(wp), INTENT(in) ::    &
      & wind_threshold,        &
      & wind_damage_scale,     &
      & cover_fract_pot,       &
      & prev_day_max_wind_10m, &
      & max_wind_10m

    REAL(wp), INTENT(inout) :: &
      & damaged_fract

    ! Locals
    REAL(wp) :: delta_time_yr

    delta_time_yr = 1._wp / 365._wp

    IF (      cover_fract_pot > EPSILON(1._wp)                      &
        .AND. prev_day_max_wind_10m > max_wind_10m * wind_threshold &
       ) THEN
      damaged_fract = MIN( 1._wp - EPSILON(1._wp) / cover_fract_pot,                                          &
        &                 delta_time_yr * wind_damage_scale * prev_day_max_wind_10m ** 3._wp / max_wind_10m   &
        &               )
    ENDIF

  END SUBROUTINE broken_woody_fract_jsbach

  FUNCTION get_relative_humidity_air(nc, air_moisture, air_temperature, air_press) RESULT(relative_humidity_air)

    USE mo_physical_constants_iface, ONLY: rd  ! gas constant for dry air [J/(K*kg)]
    USE mo_physical_constants_iface, ONLY: rv  ! gas constant for water vapor [J/(K*kg)]

    INTEGER, INTENT(in)  :: nc
    REAL(wp), INTENT(in) :: &
      air_moisture(:),      & ! Specific humidity of air [kg kg-1]
      air_temperature(:),   & ! Temperature of air [K]
      air_press(:)            ! Pressure of air [Pa]

    REAL(wp) :: air_qsat(SIZE(air_moisture)), &
                relative_humidity_air(SIZE(air_moisture))

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_relative_humidity_air'

    CALL sat_specific_humidity( &
      & nc,                     & ! in
      & air_temperature,        & ! in
      & air_press,              & ! in
      & air_qsat                & ! out
      & )

        IF (ANY( air_qsat(:) > 0.99_wp*HUGE(1._wp) )) THEN
           CALL finish(TRIM(routine), 'lookup table overflow.')
        END IF

    relative_humidity_air(:) = 100._wp * air_moisture(:)                                       &
                               * ( (1._wp - rd / rv)  * air_qsat(:)  + (rd / rv) )             &
                               / ( ((1._wp - rd / rv) * air_moisture(:) + (rd / rv)) * air_qsat(:) )

  END FUNCTION get_relative_humidity_air

  ! R: TBD: I copied this subroutine from MO_JSB4_FORCING, because
  !    this mo is not compiled with icon/echam.
  !    Reiner has to decide if we change config or the place of this procedure.
  !    I made some changes because the dimensions (vector lenghts) were unnecessary
  !    cryptic. We could even skip "nc" if we would change the DO loop below...
  !>
  !! Returns saturation specific humidity for given temperature and pressure
  !!
  SUBROUTINE sat_specific_humidity(nc, temp, pressure, qsat)

    USE mo_phy_schemes, ONLY: qsat_water

    ! Returns saturation specific humidity for given temperature and pressure (at some atmospheric level or at surface)
    ! Uses Eq. 2.27 of Fundamentals of Atmospheric Modeling for the saturated case, but the saturation vapor pressure
    ! of water over a liquid surface resp. ice (see pp 32-34 of Ref.) is computed as in ECHAM5 (mo_convect_tables)

    INTEGER,  INTENT(in) :: nc
    REAL(wp), INTENT(IN) :: temp(:)      ! Air temperature at level [K]
    REAL(wp), INTENT(IN) :: pressure(:)  ! Pressure at level [Pa]
    REAL(wp), INTENT(OUT):: qsat(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':sat_specific_humidity'

    INTEGER :: i

    DO i = 1,nc
      IF (temp(i) < 50._wp .OR. temp(i) > 400._wp) THEN
        WRITE (message_text,*) &
          & 'Temperature out of bound: ', temp(i), 'K. One thing to check: consistency of land-sea mask and forcing.'
        CALL finish(TRIM(routine), TRIM(message_text))
      ELSE
        qsat(i) = qsat_water(temp(i), pressure(i), use_convect_tables=.TRUE.)
      ENDIF
    ENDDO

  END SUBROUTINE sat_specific_humidity

#endif
END MODULE mo_disturb_process
