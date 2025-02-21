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

PROGRAM test_divide_cell_mpi

  USE mo_io_units, ONLY: nerr

  IMPLICIT NONE

  WRITE (nerr, '(a)') 'MPI test skipped in no-MPI configuration.'
END PROGRAM test_divide_cell_mpi
