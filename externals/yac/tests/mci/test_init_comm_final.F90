! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

program test_init_comm_finalize

  use utest
  use yac
  use mpi

  implicit none

  integer  :: comp_id
  integer  :: npes, mype, couple_npes, comp_npes
  integer  :: ierror
  integer  :: couple_communicator, comp_communicator

  call start_test("yac_finit and yac_ffinalize")

  call MPI_INIT(ierror)

  call MPI_COMM_RANK(MPI_COMM_WORLD, mype, ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, npes, ierror)
  call test(npes == 6)

  if (mype < 5) then

    call MPI_COMM_SPLIT(MPI_COMM_WORLD, 0, 0, couple_communicator, ierror)

    call MPI_COMM_SIZE(couple_communicator, couple_npes, ierror)
    call test(couple_npes == 5)

    call yac_finit_comm( couple_communicator)

    call yac_fdef_comp(MERGE('comp_a','comp_b', mype < 3), comp_id)

    call yac_fget_comp_comm(comp_id, comp_communicator)
    call MPI_COMM_SIZE(comp_communicator, comp_npes, ierror)
    call test(comp_npes == MERGE(3, 2, mype < 3))
    call MPI_COMM_FREE(comp_communicator, ierror)

    call yac_fenddef()

    call yac_ffinalize();

  else
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, 1, 0, couple_communicator, ierror)

    call MPI_COMM_SIZE(couple_communicator, couple_npes, ierror)
    call test(couple_npes == 1)
  endif

  call MPI_COMM_FREE(couple_communicator, ierror)
  call MPI_FINALIZE(ierror)
  call stop_test
  call exit_tests

end program test_init_comm_finalize

