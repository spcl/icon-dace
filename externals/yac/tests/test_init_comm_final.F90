!!
!! @file test_init_comm_final.F90
!!
!! @copyright Copyright  (C)  2021 DKRZ, MPI-M
!!
!! @author Moritz Hanke <hanke@dkrz.de>
!!
!
! Keywords:
! Maintainer: Moritz Hanke <hanke@dkrz.de>
!             Rene Redler <rene.redler@mpimet.mpg.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yac/
!
! This file is part of YAC.
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!


#include "test_macros.inc"

program test_init_comm_finalize

  use utest
  use mo_yac_finterface
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

