!> Contains the routines for the fuel processes
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
MODULE mo_fuel_process
#ifndef __NO_JSBACH__

  USE mo_kind,      ONLY: wp

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: calc_fuel_jsbach

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_process'

CONTAINS

  ! Calculate fuel (JSBACH3 algorithm)
  SUBROUTINE calc_fuel_jsbach( &
    & c_acid_ag1_ta,                          & ! in
    & c_water_ag1_ta,                         & ! in
    & c_ethanol_ag1_ta,                       & ! in
    & c_nonsoluble_ag1_ta,                    & ! in
    & c_acid_ag2_ta,                          & ! in
    & c_water_ag2_ta,                         & ! in
    & c_ethanol_ag2_ta,                       & ! in
    & c_nonsoluble_ag2_ta,                    & ! in
    & fuel                                    & ! out
    & )

    ! Input Arguments
    REAL(wp), INTENT(in)  :: c_acid_ag1_ta(:)
    REAL(wp), INTENT(in)  :: c_water_ag1_ta(:)
    REAL(wp), INTENT(in)  :: c_ethanol_ag1_ta(:)
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag1_ta(:)
    REAL(wp), INTENT(in)  :: c_acid_ag2_ta(:)
    REAL(wp), INTENT(in)  :: c_water_ag2_ta(:)
    REAL(wp), INTENT(in)  :: c_ethanol_ag2_ta(:)
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag2_ta(:)

    ! Output Arguments
    REAL(wp), INTENT(OUT)   :: fuel(:)

    ! Local Variables
    INTEGER :: ic, nc

    nc = SIZE(c_acid_ag1_ta)

    !$ACC PARALLEL LOOP DEFAULT(PRESENT) GANG VECTOR ASYNC(1)
    DO ic = 1, nc
      fuel(ic) =  c_acid_ag1_ta(ic)    + c_water_ag1_ta(ic)       &
        &       + c_ethanol_ag1_ta(ic) + c_nonsoluble_ag1_ta(ic)  &
        &       + c_acid_ag2_ta(ic)    + c_water_ag2_ta(ic)       &
        &       + c_ethanol_ag2_ta(ic) + c_nonsoluble_ag2_ta(ic)
    END DO
    !$ACC END PARALLEL LOOP

  END SUBROUTINE calc_fuel_jsbach

#endif
END MODULE mo_fuel_process
