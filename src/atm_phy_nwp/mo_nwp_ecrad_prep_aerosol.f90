!
! This module prepares aerosol climatologies in a format that can be used by ecRad
!
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

!----------------------------
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

!----------------------------

MODULE mo_nwp_ecrad_prep_aerosol

  USE mo_kind,                   ONLY: wp
  USE mo_exception,              ONLY: finish
  USE mo_fortran_tools,          ONLY: assert_acc_host_only, assert_acc_device_only,t_ptr_2d
  USE mo_impl_constants,         ONLY: n_camsaermr

  IMPLICIT NONE

  PRIVATE

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_nwp_ecrad_prep_aerosol'

END MODULE mo_nwp_ecrad_prep_aerosol
