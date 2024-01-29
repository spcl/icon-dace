!> Contains the memory class for the fuel process.
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
MODULE mo_fuel_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_fuel_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 5

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_fuel_memory

   TYPE(t_jsb_var_real2d) ::   &
      fuel                       ! R: Normalerweise wäre dies eine lokale Variable, weil sie nur für das veg Tile berechnet wird
                                 ! und nur für Feuer benötigt wird. Aber sie wird 1. nach dem Zeitschritt wieder verwendet und
                                 ! sollte 2. im Output sein, da sie nicht nur technisch wichtig ist!
  CONTAINS
    PROCEDURE :: Init => Init_fuel_memory
  END TYPE t_fuel_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_fuel_memory_class'

CONTAINS

  SUBROUTINE Init_fuel_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_fuel_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),     INTENT(in)    :: prefix
    CHARACTER(len=*),     INTENT(in)    :: suffix
    INTEGER,              INTENT(in)    :: lct_ids(:)
    INTEGER,              INTENT(in)    :: model_id

    TYPE(t_jsb_grid),  POINTER :: hgrid          ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface        ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_fuel_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL mem%Add_var( 'fuel', mem%fuel,                                                                      &
      & hgrid, surface,                                                                                           &
      & t_cf('fuel', 'mol(C)/m^2 (potential vegetation)', 'Amount of carbon fuel for burned area calculations.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                                        &
      & prefix, suffix,                                                                                           &
      & output_level=BASIC,                                                                                       &
      & lrestart=.FALSE., initval_r=0.0_wp ) ! initval taken from JSBACH3

  END SUBROUTINE Init_fuel_memory

#endif
END MODULE mo_fuel_memory_class
