!> Contains the memory class for the <PROCESS_NAME_LOWER_CASE> process.
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
! ===============================================================================================================================
! === THIS IS A TEMPLATE. PLEASE REPLACE ...                                                                                  ===
! === ... <PROCESS_NAME_LOWER_CASE> with your process name e.g. "seb" for surface energy balance                              ===
! ===                                                                                                                         ===
! === Then go through the code line by line to adapt it to your needs:                                                        ===
! === !X marks lines with examples (e.g:), templates (if necessary:) and implementations points (Implementation:) that you    ===
! === have to adapt to your process.                                                                                          ===
! ===                                                                                                                         ===
! === Finally delete this header section.                                                                                     ===
! ===============================================================================================================================
!>#### Contains memory definitions for <PROCESS_NAME_LOWER_CASE>
!>
MODULE mo_<PROCESS_NAME_LOWER_CASE>_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_util, ONLY: One_of

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: !X e.g: LAND_TYPE, LAKE_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_varlist,            ONLY: BASIC, MEDIUM, FULL
  USE mo_jsb_physical_constants, ONLY: !X e.g: tmelt

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_<PROCESS_NAME_LOWER_CASE>_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_<PROCESS_NAME_LOWER_CASE>_memory

    ! Common variables independent of lct_type
    TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & t,                    & !< Surface temperature    [K]


    ! Additional variables for LAND lct_type
    TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & s_star,               & !< Surface dry static energy (s^star, see manual)    [m2 s-2]


    ! Additional variables for LAKE lct_type
    TYPE(t_jsb_var_real2d) :: &
      !X add your 2D variables here e.g: & t_wtr,                & !< Lake surface temperature (water)     [K]

  CONTAINS
    PROCEDURE :: Init => Init_<PROCESS_NAME_LOWER_CASE>_memory
  END TYPE t_<PROCESS_NAME_LOWER_CASE>_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_<PROCESS_NAME_LOWER_CASE>_memory_class'

CONTAINS

  SUBROUTINE Init_<PROCESS_NAME_LOWER_CASE>_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, TSTEP_CONSTANT, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_<PROCESS_NAME_LOWER_CASE>_memory), INTENT(inout) :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid

    INTEGER :: table
    TYPE(t_grib2) :: grib2_desc

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_<PROCESS_NAME_LOWER_CASE>_memory'

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    !X add your common variables here e.g:
    !X  CALL mem%Add_var( 'sfc_temp', mem%sfc_temp,                                      &
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp', 'K', 'surface temperature'),                                &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=280.0_wp )

    !X  CALL mem%Add_var( 't_old', mem%t_old,                                            &
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp_old', 'K', 'surface temperature at previous timestep'),       &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=280.0_wp )

    ! Additional variables for LAND lct
    !X add your LAND-lct variables here e.g:
    !X IF (One_of(LAND_TYPE, lct_ids(:)) > 0) THEN

    !X  CALL mem%Add_var( 's_sfc', mem%s_star,                                                 &
    !X    & hgrid, surface,                                                                    &
    !X    & t_cf('sfc_dry_static_energy', 'm2 s-2', 'surface dry static energy'),              &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
    !X    & prefix, suffix,                                                                    &
    !X    & initval_r=2.9E5_wp )

    !X  CALL mem%Add_var( 'qsat_sfc', mem%qsat_star,                                           &
    !X    & hgrid, surface,                                                                    &
    !X    & t_cf('sfc_qsat', 'm2 s-2', 'surface specific humidity at saturation'),             &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                 &
    !X    & prefix, suffix,                                                                    &
    !X    & initval_r=0.0075_wp )

    !X END IF

    ! Additional variables for lakes
    !X add your LAKE-lct variables here e.g:
    !XIF (One_of(LAKE_TYPE, lct_ids(:)) > 0) THEN

    !X  CALL mem%Add_var( 't_wtr', mem%t_wtr,                                            &
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp_wtr', 'K', 'lake surface temperature wtr '),                  &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=280.0_wp )

    !X  CALL mem%Add_var( 't_ice', mem%t_ice,                                            &
    !X    & hgrid, surface,                                                              &
    !X    & t_cf('sfc_temp_ice', 'K', 'lake surface temperature ice '),                  &
    !X    & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),           &
    !X    & prefix, suffix,                                                              &
    !X    & initval_r=tmelt )

    !X END IF

  END SUBROUTINE Init_<PROCESS_NAME_LOWER_CASE>_memory

#endif
END MODULE mo_<PROCESS_NAME_LOWER_CASE>_memory_class
