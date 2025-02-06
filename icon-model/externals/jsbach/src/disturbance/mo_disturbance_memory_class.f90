!> Contains the memory class for the disturb process.
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
MODULE mo_disturb_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d

  dsl4jsb_Use_processes WLCC_, FLCC_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_disturb_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 40

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_disturb_memory

  TYPE(t_jsb_var_real2d) ::   &
    cconservation_fire,       &
    cconservation_wind,       &
    burned_fract,             & ! burned fraction of tile
    damaged_fract,            & ! windbreak fraction of tile
    q_rel_air_climbuf,        & ! Relative humidity of lowest air layer smoothed in time.
    q_rel_air_climbuf_yday,   & ! Relative humidity of lowest air layer smoothed in time at the end of the last day.
                                ! This is only needed for cbalone forcing!
    prev_day_max_wind_10m,    & ! for windbreak
    max_wind_10m,             & ! for windbreak
    max_wind_10m_act            ! for windbreak

  CONTAINS
    PROCEDURE :: Init => Init_disturb_memory
  END TYPE t_disturb_memory


  CHARACTER(len=*), PARAMETER :: modname = 'mo_disturb_memory_class'

CONTAINS

  SUBROUTINE Init_disturb_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: BASIC !, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    CLASS(t_disturb_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id

    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid           ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface         ! Vertical grids
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_disturb_memory'

    IF (lct_ids(1) > 0)  CONTINUE ! avoid compiler warning about dummy argument not being used
    IF (model_id > 0)    CONTINUE ! avoid compiler warning about dummy argument not being used

    model => Get_model(model_id)
    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    CALL mem%Add_var( 'burned_fract', mem%burned_fract,                                 &
      & hgrid, surface,                                                                 &
      & t_cf('burned_fract', 'frac', 'Fraction of tile that is burned.'),               &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),              &
      & prefix, suffix,                                                                 &
      & output_level=BASIC,                                                             &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'damaged_fract', mem%damaged_fract,                                &
      & hgrid, surface,                                                                  &
      & t_cf('damaged_fract', 'frac', 'Fraction of tile that is damaged by windbreak.'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),               &
      & prefix, suffix,                                                                  &
      & output_level=BASIC,                                                              &
      & loutput = .TRUE., lrestart=.TRUE., initval_r=0.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'q_air_clim', mem%q_rel_air_climbuf,                                            &
      & hgrid, surface,                                                                               &
      & t_cf('q_rel_air_climbuf', '-', 'Relative humiditiy of lowest air layer smoothed in time.'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                            &
      & prefix, suffix,                                                                               &
      & lrestart=.TRUE., initval_r=50.0_wp ) ! initval taken from JS3; lrestart=true only necessary if we restart within a day

    CALL mem%Add_var( 'q_air_clim_yday', mem%q_rel_air_climbuf_yday,                                    &
      & hgrid, surface,                                                                                 &
      & t_cf('q_rel_air_climbuf_yday', '-', 'Relative humiditiy of lowest air layer smoothed in time    &
                                               & at the end of the last day.'),                         &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                              &
      & prefix, suffix,                                                                                 &
      & output_level=BASIC,                                                                             &
      & lrestart=.FALSE., initval_r=50.0_wp ) ! initval taken from JSBACH3

    CALL mem%Add_var( 'prev_day_max_wind_10m', mem%prev_day_max_wind_10m,         &
      & hgrid, surface,                                                           &
      & t_cf('prev_day_max_wind_10m', '-', 'Previous day maximum wind speed.'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),        &
      & prefix, suffix,                                                           &
      & lrestart=.FALSE., initval_r=0.0_wp )

    CALL mem%Add_var( 'max_wind_10m', mem%max_wind_10m,                                                &
      & hgrid, surface,                                                                                &
      & t_cf('max_wind_10m', '-', 'Maximum wind speed smoothed in time at the end of the last day..'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                             &
      & prefix, suffix,                                                                                &
      & lrestart=.TRUE., initval_r=0.0_wp )

    CALL mem%Add_var( 'max_wind_10m_act', mem%max_wind_10m_act,             &
      & hgrid, surface,                                                     &
      & t_cf('max_wind_10m_act', '-', 'This days maximum wind speed.'),     &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
      & prefix, suffix,                                                     &
      & lrestart=.TRUE., initval_r=0.0_wp )

    IF (.NOT. model%processes(WLCC_)%p%config%active) THEN
      CALL mem%Add_var( 'cconservation_wind', mem%cconservation_wind,         &
        & hgrid, surface,                                                     &
        & t_cf('cconservation_wind', 'mol(C) m-2(canopy)',                    &
        &   'Test for c conservation on relocation for wind damaged area.'),  &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & loutput = .TRUE., lrestart=.FALSE., initval_r=0.0_wp )
    END IF

    IF (.NOT. model%processes(FLCC_)%p%config%active) THEN
      CALL mem%Add_var( 'cconservation_fire', mem%cconservation_fire,         &
        & hgrid, surface,                                                     &
        & t_cf('cconservation_fire', 'mol(C) m-2(canopy)',                    &
        &  'Test for c conservation on relocation for burned area.'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),  &
        & prefix, suffix,                                                     &
        & loutput = .TRUE., lrestart=.FALSE., initval_r=0.0_wp )
    END IF

  END SUBROUTINE Init_disturb_memory

#endif
END MODULE mo_disturb_memory_class
