! Contains the implementation of interpolation and reconstruction
! routines used by the shallow water model, including the RBF
! reconstruction routines.
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

MODULE mo_intp
  !-------------------------------------------------------------------------
  !
  USE mo_icon_interpolation_vector
  USE mo_icon_interpolation_scalar

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: verts2edges_scalar
  PUBLIC :: cells2edges_scalar
  PUBLIC :: edges2verts_scalar
  PUBLIC :: edges2cells_scalar
  PUBLIC :: edges2cells_vector
  PUBLIC :: cells2verts_scalar
  PUBLIC :: verts2cells_scalar
  PUBLIC :: cell_avg

END MODULE mo_intp
