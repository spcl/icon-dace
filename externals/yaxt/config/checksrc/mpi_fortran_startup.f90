!>
!! @file mpi_fortran_startup.f90
!! @brief Test startup of MPI-parallel Fortran programs
!!
!! @copyright Copyright  (C)  2023 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
! acx_mpirun_num_tasks=2
PROGRAM mpi_hello
!  INCLUDE "mpif.h"
  USE mpi
  USE iso_c_binding

  INTERFACE
    SUBROUTINE posix_exit(code) BIND(c, name="exit")
      IMPORT :: c_int
      INTEGER(c_int), VALUE, INTENT(in) :: code
    END SUBROUTINE posix_exit
  END INTERFACE

  INTEGER :: ierror, comm_size, comm_rank
  INTEGER :: procnum, acx_mpirun_num_tasks
  CALL mpi_init(ierror)
  CALL handle_mpi_error(ierror)
  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  CALL handle_mpi_error(ierror)
  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  CALL handle_mpi_error(ierror)
  acx_mpirun_num_tasks = 2
  procnum = 1
  CALL mpi_allreduce(mpi_in_place, procnum, 1, mpi_integer, &
       mpi_sum, mpi_comm_world, ierror)
  IF (procnum /= acx_mpirun_num_tasks) THEN
    CALL mpi_abort(mpi_comm_world, 1, ierror)
  END IF
  CALL mpi_finalize(ierror)
  CALL handle_mpi_error(ierror)
CONTAINS
  SUBROUTINE handle_mpi_error(ierror)
    INTEGER, INTENT(in) :: ierror

    INTEGER :: msg_len, msg_ierror
    CHARACTER(len=MPI_MAX_ERROR_STRING) :: msg

    IF (ierror /= MPI_SUCCESS) THEN
      CALL mpi_error_string(ierror, msg, msg_len, msg_ierror)
      WRITE (0, *) msg(1:msg_len)
      FLUSH(0)
      CALL posix_exit(1_c_int)
    END IF
  END SUBROUTINE handle_mpi_error

END PROGRAM mpi_hello
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
