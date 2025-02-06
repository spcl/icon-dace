! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

program test_init_finalize

  use utest
  use yac
  use mpi
  implicit none

  character(len=YAC_MAX_CHARLEN) :: comp_name

  integer  :: comp_id
  integer  :: npes, mype
  integer  :: ierror
  integer  :: local_communicator

  call start_test("yac_finit and yac_ffinalize")

  call yac_finit ( )

  comp_name = 'ICON-ocean'
  call yac_fdef_comp ( comp_name, comp_id )

  call yac_fget_comp_comm ( comp_id, local_communicator )

  call MPI_COMM_SIZE ( local_communicator, npes, ierror )
  call MPI_COMM_RANK ( local_communicator, mype, ierror )
  print *, mype, ':', ' we are running on ', npes, ' processors.'

  call yac_ffinalize ( )

  call stop_test

  call exit_tests

end program test_init_finalize

