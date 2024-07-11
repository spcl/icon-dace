!> pplcc (pre- and post lcc) memory class
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
!>#### Could contain memory definitions for pplcc
!>
!> currently of no use beyond infrastructural requirements
!>
MODULE mo_pplcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_kind, ONLY: wp
  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real1d, t_jsb_var_real2d
  USE mo_jsb_varlist,            ONLY: BASIC !, MEDIUM, FULL

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_pplcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_pplcc_memory
    TYPE(t_jsb_var_real2d) :: &
      tree_fract,             & ! land fraction covered with tree PFTs
      shrub_fract,            & ! land fraction covered with shrub PFTs
      grass_fract,            & ! land fraction covered with grass PFTs
      crop_fract,             & ! land fraction covered with crop PFTs
      pasture_fract,          & ! land fraction covered with pasture PFTs
      C3pft_fract,            & ! land fraction covered with C3 PFTs
      C4pft_fract,            & ! land fraction covered with C4 PFTs
      baresoil_fract            ! bare land fraction

    ! Diagnostic global sums e.g. for monitoring (only available with ICON)
    TYPE(t_jsb_var_real1d) :: &
      tree_area_gsum,         & ! global land area covered with tree PFTs
      shrub_area_gsum,        & ! global land area covered with shrub PFTs
      grass_area_gsum,        & ! global land area covered with grass PFTs
      crop_area_gsum,         & ! global land area covered with crop PFTs
      pasture_area_gsum,      & ! global land area covered with pasture PFTs
      C3pft_area_gsum,        & ! global land area covered with C3 PFTs
      C4pft_area_gsum,        & ! global land area covered with C4 PFTs
      baresoil_area_gsum        ! global bare land area

  CONTAINS
    PROCEDURE :: Init => Init_pplcc_memory
  END TYPE t_pplcc_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_pplcc_memory_class'

CONTAINS

  SUBROUTINE Init_pplcc_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_pplcc_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_pplcc_memory'
    ! -------------------------------------------------------------------------------------------------- !

    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    ! Area fractions of selected groups of PFTs - only calculated on box tile
    CALL mem%Add_var('tree_fract', mem%tree_fract,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('tree_fract', '1', 'Tree fraction relative to the grid cell area'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

    CALL mem%Add_var('shrub_fract', mem%shrub_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('shrub_fract', '1', 'Shrub fraction relative to the grid cell area'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

    CALL mem%Add_var('grass_fract', mem%grass_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('grass_fract', '1', 'Grass fraction relative to the grid cell area'),        &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

    CALL mem%Add_var('crop_fract', mem%crop_fract,                                        &
      & hgrid, surface,                                                                   &
      & t_cf('crop_fract', '1', 'Crop fraction relative to the grid cell area'),          &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

    CALL mem%Add_var('pasture_fract', mem%pasture_fract,                                  &
      & hgrid, surface,                                                                   &
      & t_cf('pasture_fract', '1', 'Pasture fraction relative to the grid cell area'),    &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

    CALL mem%Add_var('C3pft_fract', mem%C3pft_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('C3pft_fract', '1', 'Fraction of C3-PFTs relative to the grid cell area'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

    CALL mem%Add_var('C4pft_fract', mem%C4pft_fract,                                      &
      & hgrid, surface,                                                                   &
      & t_cf('C4pft_fract', '1', 'Fraction of C4-PFTs relative to the grid cell area'),   &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

    CALL mem%Add_var('baresoil_fract', mem%baresoil_fract,                                &
      & hgrid, surface,                                                                   &
      & t_cf('baresoil_fract', '1', 'Bare soil fraction relative to the grid cell area'), &
      & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),                &
      & prefix, suffix,                                                                   &
      & loutput=.FALSE., output_level=BASIC, initval_r=0.0_wp)

#ifdef __ICON__
    ! Diagnostic global areas for experiment monitoring
    !       (1d stream variables are not supported with echam)
    IF ( TRIM(suffix) == 'box' ) THEN
      CALL mem%Add_var('tree_area_gsum', mem%tree_area_gsum,                     &
        & hgrid, surface,                                                        &
        & t_cf('tree_area_gsum', 'Mio km2', 'Global tree area'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
      CALL mem%Add_var('shrub_area_gsum', mem%shrub_area_gsum,                   &
        & hgrid, surface,                                                        &
        & t_cf('shrub_area_gsum', 'Mio km2', 'Global shrub area'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
      CALL mem%Add_var('grass_area_gsum', mem%grass_area_gsum,                   &
        & hgrid, surface,                                                        &
        & t_cf('grass_area_gsum', 'Mio km2', 'Global grass area'),               &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
      CALL mem%Add_var('crop_area_gsum', mem%crop_area_gsum,                     &
        & hgrid, surface,                                                        &
        & t_cf('crop_area_gsum', 'Mio km2', 'Global crop area'),                 &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
      CALL mem%Add_var('pasture_area_gsum', mem%pasture_area_gsum,               &
        & hgrid, surface,                                                        &
        & t_cf('pasture_area_gsum', 'Mio km2', 'Global pasture area'),           &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
      CALL mem%Add_var('baresoil_area_gsum', mem%baresoil_area_gsum,             &
        & hgrid, surface,                                                        &
        & t_cf('baresoil_area_gsum', 'Mio km2', 'Global baresoil area'),         &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
      CALL mem%Add_var('C3pft_area_gsum', mem%C3pft_area_gsum,                   &
        & hgrid, surface,                                                        &
        & t_cf('C3pft_area_gsum', 'Mio km2', 'Global area of C3 PFTs'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
      CALL mem%Add_var('C4pft_area_gsum', mem%C4pft_area_gsum,                   &
        & hgrid, surface,                                                        &
        & t_cf('C4pft_area_gsum', 'Mio km2', 'Global area of C4 PFTs'),          &
        & t_grib1(table, 255, grib_bits), t_grib2(255, 255, 255, grib_bits),     &
        & prefix, suffix,                                                        &
        & lrestart=.FALSE., initval_r=0.0_wp)
    END IF
#endif

  END SUBROUTINE Init_pplcc_memory

#endif
END MODULE mo_pplcc_memory_class
