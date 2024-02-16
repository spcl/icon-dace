!>
!! @file test_redirstdout.f90
!!
!! @copyright Copyright  (C)  2021 DKRZ, MPI-M
!!
!! @author Rene Redler <rene.redler@mpimet.mpg.de>
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



! ------------------------------------------------------------------------

#include "test_macros.inc"

! ------------------------------------------------------------------------

PROGRAM test_redirstdout

  USE utest

  USE mpi

  USE mo_yac_finterface, ONLY: yac_redirstdout
  IMPLICIT NONE

  INTEGER :: ierror, comm_rank, comm_size


  INTEGER :: parallel

  CALL start_test("test_redirstdout")
  CALL mpi_init (ierror)
  CALL mpi_comm_rank(MPI_COMM_WORLD, comm_rank, ierror)
  CALL mpi_comm_size(MPI_COMM_WORLD, comm_size, ierror)

  IF (comm_rank == 0) CALL delete_output_files(comm_size)
  CALL mpi_barrier(MPI_COMM_WORLD, ierror)

  parallel = get_parallel()

  CALL yac_redirstdout("test_redirstdout", parallel, comm_rank, comm_size)

  CALL mpi_barrier(MPI_COMM_WORLD, ierror)

  CALL check_output_files(parallel, comm_size)

  CALL mpi_finalize (ierror)

  IF (comm_rank == 0) CALL delete_output_files(comm_size)

  CALL stop_test
  CALL exit_tests

CONTAINS

  FUNCTION get_parallel()

    IMPLICIT NONE

    INTEGER :: get_parallel

    INTEGER :: arg_len
    INTEGER, PARAMETER :: max_opt_arg_len = 80
    CHARACTER(max_opt_arg_len) :: optarg

    CALL check_test(COMMAND_ARGUMENT_COUNT() /= 1, &
                    'wrong number of arguments (has to be 1)')
    CALL GET_COMMAND_ARGUMENT(1, optarg, arg_len)
    CALL check_test(arg_len /= 1, 'wrong argument length (has to be 1)')
    CALL check_test(optarg(1:1) /= 'T' .AND. optarg(1:1) /= 'F', &
                    'wrong argument (has to be "T" or "F")')

    get_parallel = MERGE(1, 0, optarg(1:1) == 'T')

  END FUNCTION get_parallel

  SUBROUTINE check_test (test_value, err_string)

    USE mpi, ONLY : mpi_abort, MPI_COMM_WORLD
    USE utest

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: test_value
    CHARACTER(LEN=*) :: err_string

    INTEGER :: ierror

    IF (.NOT. test_value) RETURN

    WRITE (0,*) err_string

    CALL test ( .FALSE. )
    CALL stop_test
    CALL exit_tests
    CALL mpi_abort ( MPI_COMM_WORLD, 999, ierror )

  END SUBROUTINE check_test

  SUBROUTINE delete_output_files(comm_size)

    USE, INTRINSIC :: iso_c_binding, only : c_null_char

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: comm_size

    INTEGER :: i
    CHARACTER(128) :: filename

    INTERFACE

      SUBROUTINE C_UNLINK ( path ) BIND ( c, name='unlink' )

        USE, INTRINSIC :: iso_c_binding, only : c_char

        CHARACTER(KIND=c_char), DIMENSION(*) :: path

      END SUBROUTINE C_UNLINK

    END INTERFACE

    CALL C_UNLINK(TRIM('test_redirstdout.log') // c_null_char)
    CALL C_UNLINK(TRIM('test_redirstdout.err') // c_null_char)

    DO i = 0, comm_size - 1
      WRITE (filename, "('test_redirstdout.',I1)") i
      CALL C_UNLINK(TRIM(filename) // c_null_char)
      WRITE (filename, "('test_redirstdout.err.',I1)") i
      CALL C_UNLINK(TRIM(filename) // c_null_char)
    END DO

  END SUBROUTINE delete_output_files

  SUBROUTINE check_output_files(parallel, comm_size)

    INTEGER, INTENT(IN) :: parallel, comm_size

    INTEGER :: i
    LOGICAL :: file_exists
    CHARACTER(128) :: filename

    IF (parallel == 0) THEN

      INQUIRE(FILE='test_redirstdout.log', EXIST=file_exists)
      CALL check_test(.NOT. file_exists, 'test_redirstdout.log does not exist')

      INQUIRE(FILE='test_redirstdout.err', EXIST=file_exists)
      CALL check_test(.NOT. file_exists, 'test_redirstdout.err does not exist')

      DO i = 0, comm_size - 1

        WRITE (filename, "('test_redirstdout.',I1)") i
        INQUIRE(FILE=filename, EXIST=file_exists)
        CALL check_test(file_exists, 'test_redirstdout.* should not exist')

        WRITE (filename, "('test_redirstdout.err.',I1)") i
        INQUIRE(FILE=filename, EXIST=file_exists)
        CALL check_test(file_exists, 'test_redirstdout.err.* should not exist')
      END DO

    ELSE

      INQUIRE(FILE='test_redirstdout.log', EXIST=file_exists)
      CALL check_test(file_exists, 'test_redirstdout.log should not exist')

      INQUIRE(FILE='test_redirstdout.err', EXIST=file_exists)
      CALL check_test(file_exists, 'test_redirstdout.err should not exist')

      DO i = 0, comm_size - 1

        WRITE (filename, "('test_redirstdout.',I1)") i
        INQUIRE(FILE=filename, EXIST=file_exists)
        CALL check_test(.NOT. file_exists, 'test_redirstdout.* does not exist')

        WRITE (filename, "('test_redirstdout.err.',I1)") i
        INQUIRE(FILE=filename, EXIST=file_exists)
        CALL check_test(.NOT. file_exists, 'test_redirstdout.err.* does not exist')
      END DO

    END IF

  END SUBROUTINE check_output_files

END PROGRAM test_redirstdout
