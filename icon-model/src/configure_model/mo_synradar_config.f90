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

! @brief Setup for synthetic radar data on the model grid
!
! configuration setup for synthetic radar data on the model grid

MODULE mo_synradar_config
  
  USE mo_io_units,                ONLY: filename_max
  USE radar_dbzcalc_params_type,  ONLY: t_dbzcalc_params

  IMPLICIT NONE
  PUBLIC

  !--------------------------------------------------------------------------
  ! Namelist parameters
  !--------------------------------------------------------------------------

  ! Meta data for reflectivity computations (DBZ, DBZ850, DBZ_CMAX, etc.) on the model grid by using advanced methods
  !  from EMVORADO (Mie-scattering, T-matrix):
  TYPE(t_dbzcalc_params)        :: synradar_meta
  CHARACTER(LEN=filename_max) :: ydir_mielookup_read
  CHARACTER(LEN=filename_max) :: ydir_mielookup_write


END MODULE mo_synradar_config
