! This module provides data structures for reading a NetCDF file in a distributed way.
!
!
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

MODULE mo_read_netcdf_types

  USE mo_communication_types, ONLY: t_comm_pattern

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_distrib_read_data

  TYPE t_distrib_read_data
    INTEGER :: basic_data_index = -1
    CLASS(t_comm_pattern), POINTER :: pat => NULL()
  END TYPE t_distrib_read_data

END MODULE mo_read_netcdf_types
