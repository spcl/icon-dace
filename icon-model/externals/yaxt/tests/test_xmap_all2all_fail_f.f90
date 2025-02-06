!>
!! @file test_xmap_all2all_fail_f.f90
!!
!! @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!
!
! Keywords:
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
#include "fc_feature_defs.inc"
PROGRAM test_xmap_all2all_fail
  USE iso_c_binding, ONLY: c_int
  USE mpi
  USE ftest_common, ONLY: init_mpi, finish_mpi, test_abort
  USE test_idxlist_utils, ONLY: test_err_count
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_int_kind, &
       xt_stripe, xt_idxlist, xt_xmap
  USE test_xmap_common, ONLY : test_self_xmap_construct
  IMPLICIT NONE
  INTERFACE
    FUNCTION xmap_new_fail3(src_idxlist, dst_idxlist, comm) RESULT(xmap)
#if __INTEL_COMPILER == 1500
      USE yaxt, ONLY: xt_xmap, xt_idxlist
#else
      IMPORT :: xt_xmap, xt_idxlist
#endif
      TYPE(xt_xmap) :: xmap
      TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
      INTEGER, INTENT(in) :: comm
    END FUNCTION xmap_new_fail3
  END INTERFACE
  INTEGER, PARAMETER :: xi = xt_int_kind
  CHARACTER(len=*), PARAMETER :: filename = 'test_xmap_all2all_fail_f.f90'
  INTEGER :: my_rank, ierror, list_size
  CALL init_mpi
  CALL xt_initialize(mpi_comm_world)
  CALL mpi_comm_rank(mpi_comm_world, my_rank, ierror)
  CALL parse_options
  CALL test_xmap1(list_size, mpi_comm_world)

  IF (test_err_count() /= 0) &
       CALL test_abort("non-zero error count!", filename, __LINE__)
  CALL xt_finalize
  CALL finish_mpi
CONTAINS
  SUBROUTINE test_xmap1(num_idx, comm)
    INTEGER, INTENT(in) :: num_idx, comm
    TYPE(xt_stripe) :: src_stripe(1), dst_stripe(1)

    ! soruce index list
    src_stripe(1)%nstrides = INT(num_idx, c_int)
    src_stripe(1)%start = 1_xi + INT(my_rank, xi) * INT(num_idx, xi)
    src_stripe(1)%stride = 1_xi

    ! destination index list
    dst_stripe(1)%nstrides = INT(num_idx, c_int)
    dst_stripe(1)%start = src_stripe(1)%start + src_stripe(1)%nstrides
    dst_stripe(1)%stride = -1_xi
    ! note: this should fail because dst/src indices don't match
    CALL test_self_xmap_construct(src_stripe, dst_stripe,&
         & xmap_new_fail3, comm)
  END SUBROUTINE test_xmap1

  SUBROUTINE parse_options
    INTEGER :: i, num_cmd_args, arg_len
    INTEGER, PARAMETER :: max_opt_arg_len = 80
    CHARACTER(max_opt_arg_len) :: optarg
    num_cmd_args = COMMAND_ARGUMENT_COUNT()
    i = 1
    DO WHILE (i < num_cmd_args)
      CALL GET_COMMAND_ARGUMENT(i, optarg, arg_len)
      IF (optarg(1:2) == '-s' .AND. i < num_cmd_args .AND. arg_len == 2) THEN
        CALL GET_COMMAND_ARGUMENT(i + 1, optarg, arg_len)
        IF (arg_len > max_opt_arg_len) &
             CALL test_abort('incorrect argument to command-line option -s', &
             filename, __LINE__)
        IF (optarg(1:arg_len) == "big") THEN
          list_size = 1023
        ELSE IF (optarg(1:arg_len) == "small") THEN
          list_size = 7
        ELSE
          WRITE (0, *) 'arg to -s: ', optarg(1:arg_len)
          CALL test_abort('incorrect argument to command-line option -s', &
               filename, __LINE__)
        END IF
        i = i + 2
      ELSE
        WRITE (0, *) 'unexpected command-line argument parsing error: ', &
             optarg(1:arg_len)
        FLUSH(0)
        CALL test_abort('unexpected command-line argument', &
             filename, __LINE__)
      END IF
    END DO
  END SUBROUTINE parse_options

END PROGRAM test_xmap_all2all_fail

SUBROUTINE xfail_abort(comm, msg, source, line)
  USE iso_c_binding, ONLY: c_int
  USE mpi
  USE ftest_common, ONLY: posix_exit
  INTEGER, INTENT(in) :: comm, line
  CHARACTER(len=*), INTENT(in) :: msg, source
  INTEGER :: ierror
#ifdef XT_NEED_MPI_ABORT_WORK_AROUND
  INTEGER :: abort_msg_lun
#endif
  WRITE (0, '(4a,i0)') msg, ' at ', source, ', line ', line

#ifdef XT_NEED_MPI_ABORT_WORK_AROUND
#  if XT_NEED_MPI_ABORT_WORK_AROUND == 1
  abort_msg_lun = 0
#  elif XT_NEED_MPI_ABORT_WORK_AROUND == 2
  FLUSH(0)
  abort_msg_lun = 10
  OPEN(unit=abort_msg_lun, status='replace', &
       file='test_xmap_all2all_fail.result.txt', &
       action='write', iostat=ierror)
#  endif
  WRITE (abort_msg_lun, '(a)') 'MPI_Abort(0xdeadbeef, 3)'
  FLUSH(abort_msg_lun)
#endif

  CALL mpi_abort(comm, 3, ierror)
  CALL posix_exit(3_c_int)
END SUBROUTINE xfail_abort

FUNCTION xmap_new_fail3(src_idxlist, dst_idxlist, comm) RESULT(xmap)
  USE yaxt, ONLY: xt_xmap, xt_idxlist, xt_xmap_all2all_new, &
       xt_set_abort_handler, xt_restore_default_abort_hndl
  TYPE(xt_xmap) :: xmap
  TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
  INTEGER, INTENT(in) :: comm
  INTERFACE
    SUBROUTINE xfail_abort(comm, msg, source, line)
      INTEGER, INTENT(in) :: comm, line
      CHARACTER(len=*), INTENT(in) :: msg, source
    END SUBROUTINE xfail_abort
  END INTERFACE
  CALL xt_set_abort_handler(xfail_abort)
  xmap = xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm)
  CALL xt_restore_default_abort_hndl
END FUNCTION xmap_new_fail3

!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
