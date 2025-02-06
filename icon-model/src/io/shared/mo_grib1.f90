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

MODULE mo_grib1

  IMPLICIT NONE

  PRIVATE

  TYPE t_grib1_global
    INTEGER :: centre
    INTEGER :: subcentre
  END TYPE t_grib1_global

  TYPE t_grib1_var
    INTEGER :: table
    INTEGER :: parameter
    INTEGER :: bits
    INTEGER :: leveltype
  END type t_grib1_var

  PUBLIC :: t_grib1_global
  PUBLIC :: t_grib1_var

END MODULE mo_grib1
