# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/tests/test_divide_cell_mpi.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/tests/test_divide_cell_mpi.f90"
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

PROGRAM test_divide_cell_mpi

# 26 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/tests/test_divide_cell_mpi.f90"
  USE mo_io_units, ONLY: nerr

  IMPLICIT NONE


  WRITE (nerr, '(a)') 'MPI test skipped in no-MPI configuration.'
# 272 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/src/tests/test_divide_cell_mpi.f90"
END PROGRAM test_divide_cell_mpi
