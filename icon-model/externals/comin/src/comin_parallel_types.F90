!> @file comin_parallel_types.F90
!! @brief Data types for parallel communication.
!
!  @authors 10/2023 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_parallel_types
  USE mpi
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_comin_parallel_info

  ! max. character string length (e.g. MPI group name)
  INTEGER, PARAMETER :: MAX_GRPNAMELEN = 256

  ! data structure containing parallelization (MPI-) related data.
  TYPE :: t_comin_parallel_info
    ! MPI communicator, comprising ICON's participating PEs
    INTEGER :: host_comm  =  MPI_COMM_NULL

    ! list of MPI intra-communicators, formed by the PEs of each group
    ! (the term "group" can be identified with the set of external MPI
    ! processes that communicate with a specific plugin)
    ! plus a "root group". If the current PE is not a member of the
    ! root group, then this list contains a single entry only ("root
    ! PEs" + PEs of the local group).
    INTEGER,                       ALLOCATABLE :: mpi_comm_bilateral(:)
    CHARACTER(LEN=MAX_GRPNAMELEN), ALLOCATABLE :: component_name_bilateral(:)
  END TYPE t_comin_parallel_info

END MODULE comin_parallel_types
