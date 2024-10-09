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

! Some CDI-specific constants

MODULE mo_cdi_constants

  PUBLIC

  !------------------------------------------------!
  !  CDI constants for horizontal grid
  !------------------------------------------------!

  INTEGER, PARAMETER :: GRID_CELL   = 1
  INTEGER, PARAMETER :: GRID_VERTEX = 2
  INTEGER, PARAMETER :: GRID_EDGE   = 3

  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_CELL = 1
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_VERT = 2
  INTEGER, PARAMETER :: GRID_UNSTRUCTURED_EDGE = 3

  INTEGER, PARAMETER :: GRID_ZONAL = 16

  INTEGER, PARAMETER :: GRID_REGULAR_LONLAT    = 45

  CHARACTER(*), PARAMETER :: LONLAT_PREFIX = "latlon_internal_"
END MODULE mo_cdi_constants
