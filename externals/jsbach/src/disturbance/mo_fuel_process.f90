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
  PURE ELEMENTAL SUBROUTINE calc_fuel_jsbach( &
    & c_acid_ag1_ta,                          & ! in
    & c_water_ag1_ta,                         & ! in
    & c_ethanol_ag1_ta,                       & ! in
    & c_nonsoluble_ag1_ta,                    & ! in
    & c_acid_ag2_ta,                          & ! in
    & c_water_ag2_ta,                         & ! in
    & c_ethanol_ag2_ta,                       & ! in
    & c_nonsoluble_ag2_ta,                    & ! in
    & fract_fpc_max,                          & ! in
    & fuel                                    & ! out
    & )

    ! Input Arguments
    REAL(wp), INTENT(in)  :: c_acid_ag1_ta
    REAL(wp), INTENT(in)  :: c_water_ag1_ta
    REAL(wp), INTENT(in)  :: c_ethanol_ag1_ta
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag1_ta
    REAL(wp), INTENT(in)  :: c_acid_ag2_ta
    REAL(wp), INTENT(in)  :: c_water_ag2_ta
    REAL(wp), INTENT(in)  :: c_ethanol_ag2_ta
    REAL(wp), INTENT(in)  :: c_nonsoluble_ag2_ta
    REAL(wp), INTENT(in)  :: fract_fpc_max

    ! Output Arguments
    REAL(wp), INTENT(OUT)   :: fuel

    ! Locals

    ! ---------------------------
    ! Go

    fuel =  (   c_acid_ag1_ta    + c_water_ag1_ta       &
      &       + c_ethanol_ag1_ta + c_nonsoluble_ag1_ta  &
      &       + c_acid_ag2_ta    + c_water_ag2_ta       &
      &       + c_ethanol_ag2_ta + c_nonsoluble_ag2_ta  &
      &     ) * fract_fpc_max  ! According to fire_litter_threshold, fuel should be relative to potential vegetation area of the veg tile.
                               ! R: Note, we do not need veg_fract_correction. It was already considered before aggregating the C-pools_ta
                               !    from pft areas to this veg area. Anyhow it would not be useable on the veg tile here.

  END SUBROUTINE calc_fuel_jsbach

#endif
END MODULE mo_fuel_process
