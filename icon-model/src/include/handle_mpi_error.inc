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

#define HANDLE_MPI_ERROR_USE use mo_mpi, only: handle_mpi_error, mpi_success
#define HANDLE_MPI_ERROR(ierror, mpi_routine) IF (ierror /= mpi_success) \
call handle_mpi_error(ierror, routine, __LINE__, mpi_routine)
