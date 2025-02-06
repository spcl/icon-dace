! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

! Module containing type defintions and constructors for GRIB2 tile output
!
! GRIB2 tile output can be written by using either the official WMO tile
! templates 55, 59, ..., or the local DWD templates 40455, 40456. These
! templates do not differ only in terms of the template number but also
! wrt the naming of tile-specific keys.
! The aim of the templates below is to enable ICON to make use of either
! WMO or DWD templates for writing GRIB2 tile output.

MODULE mo_grib2_tile

  IMPLICIT NONE

  PRIVATE


  ! tile keynames
  ! depend on the actual set of tile templates in use
  TYPE t_grib2_keys_tile
    CHARACTER(len=50) :: tileClassification
    CHARACTER(len=50) :: totalNumberOfTileAttributePairs
    CHARACTER(len=50) :: numberOfUsedSpatialTiles
    CHARACTER(len=50) :: tileIndex
    CHARACTER(len=50) :: numberOfUsedTileAttributes
    CHARACTER(len=50) :: attributeOfTile
  END TYPE t_grib2_keys_tile

  ! Type for defining the set of employed tile templates 
  ! and corresponding key names.
  TYPE t_grib2_template_tile
    ! set of allowed tile templates
    ! (either offical WMO ones or ones from DWD)
    INTEGER :: tpl_inst       ! template instant, deterministic
    INTEGER :: tpl_acc        ! template acc/avg, deterministic
    INTEGER :: tpl_inst_ens   ! template instant, ensemble
    INTEGER :: tpl_acc_ens    ! template acc/avg ensemble
    !
    ! tile keynames
    TYPE(t_grib2_keys_tile):: keys
  END TYPE t_grib2_template_tile


  ! types
  PUBLIC :: t_grib2_keys_tile
  PUBLIC :: t_grib2_template_tile

  ! constructor
  PUBLIC :: grib2_keys_tile

CONTAINS

  ! constructor for derived type t_grib2_keys_tile
  FUNCTION grib2_keys_tile(str_tileClassification, str_totalNumberOfTileAttributePairs, &
    &                      str_numberOfUsedSpatialTiles, str_tileIndex,                 &
    &                      str_numberOfUsedTileAttributes, str_attributeOfTile)
    CHARACTER(len=*), INTENT(IN) :: str_tileClassification
    CHARACTER(len=*), INTENT(IN) :: str_totalNumberOfTileAttributePairs
    CHARACTER(len=*), INTENT(IN) :: str_numberOfUsedSpatialTiles
    CHARACTER(len=*), INTENT(IN) :: str_tileIndex
    CHARACTER(len=*), INTENT(IN) :: str_numberOfUsedTileAttributes
    CHARACTER(len=*), INTENT(IN) :: str_attributeOfTile
    TYPE(t_grib2_keys_tile) :: grib2_keys_tile

    grib2_keys_tile%tileClassification              = str_tileClassification
    grib2_keys_tile%totalNumberOfTileAttributePairs = str_totalNumberOfTileAttributePairs
    grib2_keys_tile%numberOfUsedSpatialTiles        = str_numberOfUsedSpatialTiles
    grib2_keys_tile%tileIndex                       = str_tileIndex
    grib2_keys_tile%numberOfUsedTileAttributes      = str_numberOfUsedTileAttributes
    grib2_keys_tile%attributeOfTile                 = str_attributeOfTile
  END FUNCTION grib2_keys_tile


END MODULE mo_grib2_tile
