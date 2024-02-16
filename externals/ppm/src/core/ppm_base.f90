!>
!! @file ppm_base.f90
!! @brief Fortran 90 interface to core functions
!!
!! @copyright Copyright  (C)  2010  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords: core wrapper
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
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
! Commentary:
!
! Look into core.c for functionality common to Fortran and C parts of PPM
! Also provides extra functionality for concise expressions.
!
! Code:
!
#include "fc_feature_defs.inc"
MODULE ppm_base
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#ifdef USE_MPI
#ifndef USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
#else
  !> communicator object to use by default
  INTEGER, PARAMETER :: mpi_comm_world = 0
#endif
  !> communicator object to use by default
  INTEGER :: ppm_default_comm
  INCLUDE 'ppmcommon.inc'
  SAVE :: /ppm_f2c_data/
  PUBLIC :: ppm_default_comm, set_default_comm
  PUBLIC :: abort_ppm, assertion
  PUBLIC :: set_abort_handler, restore_default_abort_handler
  !> this should go into a wrapper module
  PUBLIC :: mpi_comm_world
#ifdef USE_MPI
  PUBLIC :: calls_to_mpi_are_allowed
#endif
CONTAINS
  !> abort operation in library, this will call the function reference
  !! assigned to PPM_abort on the C side and substitute non-provided
  !! optional dummy arguments as needed.
  !!
  !! @param msg text to write to standard error
  !! @param source string describing source file name
  !! @param line line number of caller
  !! @param comm communicator to use in PPM_abort call, defaults
  !! to ppm_default_comm if not given
  SUBROUTINE abort_ppm(msg, source, line, comm)
    CHARACTER(len=*), INTENT(in) :: source, msg
    INTEGER, INTENT(in) :: line
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER :: comm_dummy

    ! abort operation in library, this subroutine will call the
    ! function reference assigned to PPM_abort on the C side directly.
    !
    ! @param comm communicator to use in PPM_abort call, defaults
    ! to ppm_default_comm if not given
    ! @param msg text to write to standard error
    ! @param source string describing source file name
    ! @param line line number of caller
    INTERFACE
      SUBROUTINE ppm_abort(comm, msg, source, line)
        INTEGER, INTENT(in) :: comm, line
        CHARACTER(*), INTENT(in) :: msg, source
      END SUBROUTINE ppm_abort
    END INTERFACE

    IF (PRESENT(comm)) THEN
      comm_dummy = comm
    ELSE
      comm_dummy = ppm_default_comm
    END IF
    CALL ppm_abort(comm_dummy, msg, source, line)
  END SUBROUTINE abort_ppm

  !> change value of default communicator object used in library calls
  !! parallel routines in the library will assume that all
  !! participating processes are organized in this communicator
  SUBROUTINE set_default_comm(comm)
    INTEGER, INTENT(in) :: comm

    INTERFACE
      !  change default communicator object on C side
      SUBROUTINE ppm_set_default_comm(comm)
        INTEGER, INTENT(in) :: comm
      END SUBROUTINE ppm_set_default_comm
    END INTERFACE

    ppm_default_comm = comm
    CALL ppm_set_default_comm(comm)
  END SUBROUTINE set_default_comm

  !> set routine f to use as abort function which is called on ppm_abort
  SUBROUTINE set_abort_handler(f)
    INTERFACE
      SUBROUTINE f(comm, msg, source, line)
        INTEGER, INTENT(in) :: comm, line
        CHARACTER(len=*), INTENT(in) :: msg, source
      END SUBROUTINE f
      SUBROUTINE ppm_set_abort_handler(f)
        INTERFACE
          SUBROUTINE f(comm, msg, source, line)
            INTEGER, INTENT(in) :: comm, line
            CHARACTER(len=*), INTENT(in) :: msg, source
          END SUBROUTINE f
        END INTERFACE
      END SUBROUTINE ppm_set_abort_handler
    END INTERFACE
    CALL ppm_set_abort_handler(f)
  END SUBROUTINE set_abort_handler

  SUBROUTINE restore_default_abort_handler
    INTERFACE
      SUBROUTINE ppm_restore_default_abort_hndl
      END SUBROUTINE ppm_restore_default_abort_hndl
    END INTERFACE
    CALL ppm_restore_default_abort_hndl
  END SUBROUTINE restore_default_abort_handler

  !> check invariant and call abort_ppm if false
  !!
  !! @param cond invariant to test
  !! @param source string describing source file name
  !! @param line line number of caller
  !! @param msg text to use as diagnostic message
  SUBROUTINE assertion(cond, source, line, msg)
    LOGICAL, INTENT(in) :: cond
    CHARACTER(len=*), OPTIONAL, INTENT(in) :: source, msg
    INTEGER, OPTIONAL, INTENT(in) :: line
    CHARACTER(len=255) :: source_arg
    INTEGER :: line_arg
    CHARACTER(*), PARAMETER :: default_msg = 'assertion failed'

    IF (cond) RETURN
    IF (PRESENT(source)) THEN
      source_arg = source
    ELSE
      source_arg = 'unknown'
    ENDIF
    IF (PRESENT(line)) THEN
      line_arg = line
    ELSE
      line_arg = -1
    ENDIF
    IF (PRESENT(msg)) THEN
      CALL abort_ppm(msg, TRIM(source_arg), line_arg)
    ELSE
      CALL abort_ppm(default_msg, TRIM(source_arg), line_arg)
    END IF
  END SUBROUTINE assertion

#ifdef USE_MPI
  FUNCTION calls_to_mpi_are_allowed() RESULT(p)
    LOGICAL :: init_flag, finished_flag, p
    INTEGER :: ierror(2)
    init_flag = .FALSE.
    finished_flag = .FALSE.
    CALL mpi_initialized(init_flag, ierror(1))
    IF (ierror(1) == mpi_success) &
         CALL mpi_finalized(finished_flag, ierror(2))
    p = ierror(1) == MPI_SUCCESS .AND. init_flag &
         .AND. ierror(2) == MPI_SUCCESS .AND. .NOT. finished_flag
  END FUNCTION calls_to_mpi_are_allowed
#endif

END MODULE ppm_base
BLOCK DATA
  USE ppm_base, ONLY: mpi_comm_world
  INCLUDE 'ppmcommon.inc'
  SAVE :: /ppm_f2c_data/
  INTEGER :: ppm_default_comm
  DATA ppm_default_comm /mpi_comm_world/
END BLOCK DATA

!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! End:
