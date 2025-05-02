# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_surrogate_class.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_aes/tmx/mo_surrogate_class.f90"
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

! Basic class used in turbulent mixing package (tmx)

MODULE mo_surrogate_class

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_surrogate

  TYPE, ABSTRACT :: t_surrogate
  END TYPE t_surrogate

END MODULE mo_surrogate_class
