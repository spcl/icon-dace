!!
!! @file test_mpi_handshake.F90
!!
!! @copyright Copyright  (C)  2022 DKRZ, MPI-M
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

program test_mpi_handshake

  use utest
  use mo_yac_finterface
  use yaxt
  use mpi
  implicit none

  integer :: global_rank, global_size
  integer :: ierror

  call mpi_init(ierror)

  call xt_initialize(mpi_comm_world)

  call mpi_comm_rank ( mpi_comm_world, global_rank, ierror )
  call mpi_comm_size ( mpi_comm_world, global_size, ierror )

  if (global_size /= 6) then
    write ( * , * ) "wrong number of processes (should be 6)"
    call error_exit
  endif

  call test_fmpi_handshake()

  call xt_finalize()

  call mpi_finalize(ierror)

  call stop_test
  call exit_tests

contains

  subroutine test_fmpi_handshake ()

    integer :: yac_instance_id
    logical :: is_dummy
    integer :: color

    integer :: ref_group_comm, ref_main_comm
    integer :: main_comm, group_comm
    integer :: result
    integer :: group_comms(2)
    CHARACTER(len=YAC_MAX_CHARLEN) :: group_names(2)

    yac_instance_id = -1
    is_dummy = .false.
    color = mpi_undefined
    main_comm = mpi_comm_null
    group_comm = mpi_comm_null

    select case(global_rank)
      case (0:1)
         ! first two ranks are part of group A
         group_names(1) = "main"
         group_names(2) = "group a"
         CALL yac_fmpi_handshake(MPI_COMM_WORLD, group_names, group_comms)
         main_comm = group_comms(1)
         group_comm = group_comms(2)
         color = 0
      case (2:3)
         ! next two ranks are part of group B
         group_names(1) = "group b"
         group_names(2) = "main"
         CALL yac_fmpi_handshake(MPI_COMM_WORLD, group_names, group_comms)
         group_comm = group_comms(1)
         main_comm = group_comms(2)
        color = 1
      case (4)
         ! rank four does not provide a group
         group_names(1) = "main"
         CALL yac_fmpi_handshake(MPI_COMM_WORLD, group_names(1:1), group_comms(1:1))
         main_comm = group_comms(1)
      case (5)
        ! rank five does a dummy initialisation
         CALL yac_fmpi_handshake(MPI_COMM_WORLD, group_names(2:1), group_comms(2:1))
         is_dummy = .true.
    end select

    call mpi_comm_split(mpi_comm_world, color, 0, ref_group_comm, ierror)
    call mpi_comm_split( &
      mpi_comm_world, MERGE(mpi_undefined, 0, global_rank == 5), 0, &
      ref_main_comm, ierror)

    if (ref_group_comm /= mpi_comm_null) then

      call mpi_comm_compare(ref_group_comm, group_comm, result, ierror);

      call test((result == MPI_CONGRUENT) .OR. (result == MPI_IDENT))

      call MPI_Comm_free(ref_group_comm, ierror)
      call MPI_Comm_free(group_comm, ierror)

    else
      call test(group_comm == mpi_comm_null)
    end if

    if (ref_main_comm /= mpi_comm_null) then

      call mpi_comm_compare(ref_main_comm, main_comm, result, ierror);

      call test((result == MPI_CONGRUENT) .OR. (result == MPI_IDENT))

      call MPI_Comm_free(ref_main_comm, ierror)
      call MPI_Comm_free(main_comm, ierror)

    else
      call test(main_comm == mpi_comm_null)
    end if

  end subroutine test_fmpi_handshake

  subroutine error_exit ()

    use mpi, only : mpi_abort, mpi_comm_world
    use utest

    integer :: ierror

    call test ( .false. )
    call stop_test
    call exit_tests
    call mpi_abort ( mpi_comm_world, 999, ierror )

  end subroutine error_exit

end program test_mpi_handshake

