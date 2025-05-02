# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/serialization/mo_ser_manually.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/serialization/mo_ser_manually.f90"
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

! Serialize any scalar and array component in selected global derived type
! variables and that were not initialized with add_var or add_ref.
!
! HOW TO:
! 1) Add a USE-ONLY statement for each derived type variable
! 2) Add ser_component() calls to ser_manually() for each component

MODULE mo_ser_manually

# 253 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/serialization/mo_ser_manually.f90"

END MODULE mo_ser_manually
