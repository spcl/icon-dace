# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_rte_rrtmgp/mo_cloud_optics_parameters.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/atm_phy_rte_rrtmgp/mo_cloud_optics_parameters.f90"
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

! Module to provide parameters to radiation routines and avoid circular dependencies.
!
! Remarks
!   This module contains the public parameters provided by the radiation module
!   mo_radiation.

MODULE mo_cloud_optics_parameters

  USE mo_kind,   ONLY: wp

IMPLICIT NONE

  PRIVATE

  PUBLIC :: rad_perm
  
  INTEGER :: rad_perm = 1                ! Integer for perturbing random number seeds

END MODULE mo_cloud_optics_parameters
