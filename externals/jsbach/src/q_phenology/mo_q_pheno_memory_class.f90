!> QUINCY phenology process memory
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### definition and init of (memory) variables for the phenology process
!>
MODULE mo_q_pheno_memory_class
#ifndef __NO_QUINCY__

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: message, message_text, finish
  USE mo_util,                   ONLY: One_of
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory
  USE mo_jsb_lct_class,          ONLY: VEG_TYPE, LAND_TYPE
  USE mo_jsb_var_class,          ONLY: t_jsb_var_real2d, t_jsb_var_real3d

  dsl4jsb_Use_processes Q_PHENO_
  dsl4jsb_Use_config(Q_PHENO_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_q_pheno_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 20

  ! ======================================================================================================= !
  !> Type definition for pheno memory
  !>
  !>
  TYPE, EXTENDS(t_jsb_memory) :: t_q_pheno_memory

    TYPE(t_jsb_var_real2d) :: &
      & growing_season, &           !< logical: growing season for the plant [0 = FALSE, 1 = TRUE]
      & gdd, &                      !< growing degree days above 5degC [K day]
      & nd_dormance, &              !< number of days since last growing season [# days]
      & lai_max                     !< maximum lai - from site specific (i.e. local) lai data [m2 m-2]

  CONTAINS
    PROCEDURE :: Init => Init_q_pheno_memory
  END TYPE t_q_pheno_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_q_pheno_memory_class'

CONTAINS

  ! ======================================================================================================= !
  !> initialize memory (variables) for the process: q_phenology
  !>
  !>
  SUBROUTINE Init_q_pheno_memory(mem, prefix, suffix, lct_ids, model_id)
    USE mo_jsb_model_class,   ONLY: t_jsb_model
    USE mo_jsb_class,         ONLY: Get_model
    USE mo_jsb_varlist,       ONLY: BASIC, MEDIUM, FULL
    USE mo_jsb_io,            ONLY: grib_bits, t_cf, t_grib1, t_grib2, tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_q_pheno_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),        INTENT(in)            :: prefix
    CHARACTER(len=*),        INTENT(in)            :: suffix
    INTEGER,                 INTENT(in)            :: lct_ids(:)
    INTEGER,                 INTENT(in)            :: model_id
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_grid),   POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid),  POINTER :: surface                      ! Vertical grid
    TYPE(t_jsb_vgrid),  POINTER :: vgrid_canopy_q_assimi        ! Vertical grid
    INTEGER                     :: table                        ! ...
    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_q_pheno_memory'
    ! ----------------------------------------------------------------------------------------------------- !
    table                 = tables(1)
    hgrid                 => Get_grid(mem%grid_id)
    surface               => Get_vgrid('surface')
    vgrid_canopy_q_assimi => Get_vgrid('q_canopy_layer')
    ! ----------------------------------------------------------------------------------------------------- !

    ! create memory at tiles of LAND_TYPE or VEG_TYPE
    IF ( One_of(LAND_TYPE, lct_ids(:)) > 0 .OR. One_of(VEG_TYPE,  lct_ids(:)) > 0) THEN

      CALL mem%Add_var('growing_season', mem%growing_season, &
        & hgrid, surface, &
        & t_cf('growing_season', 'logical', 'logical to identify whether the plant is in the growing season'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = BASIC, &
        & loutput = .TRUE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('gdd', mem%gdd, &
        & hgrid, surface, &
        & t_cf('gdd', 'K days', 'growing degree days above 5 degC'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('nd_dormance', mem%nd_dormance, &
        & hgrid, surface, &
        & t_cf('nd_dormance', '# days', 'number of days since last growing season'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)

      CALL mem%Add_var('lai_max', mem%lai_max, &
        & hgrid, surface, &
        & t_cf('lai_max', 'm2 m-2', 'maximum lai (site specific)'), &
        & t_grib1(table, 255, grib_bits), &
        & t_grib2(255, 255, 255, grib_bits), &
        & prefix, suffix, &
        & output_level = FULL, &
        & loutput = .FALSE., &
        & lrestart = .TRUE., &
        & initval_r = 0.0_wp)
    END IF
  END SUBROUTINE Init_q_pheno_memory

#endif
END MODULE mo_q_pheno_memory_class
