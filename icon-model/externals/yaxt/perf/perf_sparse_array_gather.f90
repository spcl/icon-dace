!>
!! @file perf_sparse_array_gather.f90
!! @brief recreates whole array from variously displaced chunks
!!
!! @copyright Copyright  (C)  2021 Jörg Behrens <behrens@dkrz.de>
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
PROGRAM perf_sparse_array_gather
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_abort, &
       xt_int_kind, xt_stripe, xt_idxlist_collection_new, &
       xt_idxlist, xt_idxfsection_new, xt_idxstripes_new, xt_idxlist_delete, &
       xt_redist, xt_redist_p2p_ext_new, xt_offset_ext, &
       xt_xmap, xt_xmap_all2all_new, xt_xmap_delete, &
       xt_redist_delete, xt_redist_s_exchange, &
       xt_idxstripes_from_idxlist_new
  USE iso_c_binding, ONLY: c_int
  USE mpi
  ! older PGI compilers do not handle generic interface correctly
#if defined __PGI && ( __PGIC__ == 15 || __PGIC__ == 14 )
  USE xt_redist_real_sp, ONLY: xt_redist_s_exchange
  USE xt_redist_real_dp, ONLY: xt_redist_s_exchange
#endif
  IMPLICIT NONE

  INTEGER :: gshape(3), nparts(2), num_parts, sshape(3), gsize, &
       section_shape(3), section_start(3), section_size, ofs, ierror
  DOUBLE PRECISION, ALLOCATABLE :: gathered_data(:, :, :), &
       scattered_data(:, :, :)
  TYPE(xt_offset_ext) :: gather_ext(1)
  TYPE(xt_offset_ext), ALLOCATABLE :: scatter_exts(:,:)
  TYPE(xt_idxlist) :: gather_list
  TYPE(xt_idxlist), ALLOCATABLE :: scatter_lists(:,:)
  TYPE(xt_idxlist) :: scatter_list, scatter_list_stripes
  TYPE(xt_xmap) :: gather_xmap
  TYPE(xt_redist) :: gather_redist
  REAL :: r
  INTEGER :: i, j, k, t
  CHARACTER(len=132) :: msg
  CHARACTER(len=*), PARAMETER :: filename = 'perf_sparse_array_gather.f90'
  LOGICAL :: check_correctness_of_exchange
  INTERFACE
    SUBROUTINE posix_exit(code) BIND(c, name="exit")
      IMPORT :: c_int
      INTEGER(c_int), VALUE, INTENT(in) :: code
    END SUBROUTINE posix_exit
  END INTERFACE

  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) THEN
    WRITE (0, *) "MPI initialization failed"
    CALL POSIX_EXIT(1_c_int)
  END IF
  CALL xt_initialize(mpi_comm_world)

  check_correctness_of_exchange = .FALSE.
  CALL parse_options

  gshape(1) = 192
  gshape(2) = 96
  gshape(3) = 49
  gsize = PRODUCT(gshape)

  nparts(1) = 8
  nparts(2) = 8
  num_parts = PRODUCT(nparts)
  sshape(1) = (gsize + num_parts - 1) / num_parts * 20
  sshape(2) = nparts(1)
  sshape(3) = nparts(2)

  ALLOCATE(gathered_data(gshape(1), gshape(2), gshape(3)), &
       scattered_data(sshape(1), sshape(2), sshape(3)), &
       scatter_exts(nparts(1), nparts(2)), scatter_lists(nparts(1), nparts(2)))

  gather_ext(1) = xt_offset_ext(0, gsize, 1)

  gathered_data = -HUGE(1.0d0)

  ! perform 2d deco of cube
  ! divide gshape(1:2) by nparts(1:2)
  section_start(3) = 1
  section_shape(3) = gshape(3)
  DO t = 1, 100
    gather_list = xt_idxstripes_new(xt_stripe(0, 1, gsize))
    CALL RANDOM_NUMBER(r)
    DO j = 1, nparts(2)
      section_start(2) = (gshape(2) * (j - 1))/nparts(2) + 1
      section_shape(2) = (gshape(2) * j)/nparts(2) - section_start(2) + 1
      DO i = 1, nparts(1)
        section_start(1) = (gshape(1) * (i - 1))/nparts(1) + 1
        section_shape(1) = (gshape(1) * i)/nparts(1) - section_start(1) + 1
        section_size = PRODUCT(section_shape(:))
        ofs = INT(r * REAL(sshape(1) - 1 - section_size))
        scatter_exts(i, j) &
             = xt_offset_ext(ofs + ((i-1) + (j-1) * sshape(2)) * sshape(1), &
             section_size, 1)
        scatter_lists(i, j) = xt_idxfsection_new(0_xt_int_kind, &
             INT(gshape, xt_int_kind), &
             section_shape, INT(section_start, xt_int_kind))
        scattered_data(1:ofs, i, j) = -9.0d99
        scattered_data(ofs+1:ofs + section_size, i, j) &
             = DBLE(i + (j - 1) * nparts(1))
        scattered_data(ofs + section_size + 1:, i, j) = -9.0d99
      END DO
    END DO
    scatter_list = xt_idxlist_collection_new(scatter_lists)
    CALL xt_idxlist_delete(scatter_lists)
    scatter_list_stripes = xt_idxstripes_from_idxlist_new(scatter_list)
    CALL xt_idxlist_delete(scatter_list)
    gather_xmap = xt_xmap_all2all_new(scatter_list_stripes, gather_list, &
         mpi_comm_self)
    CALL xt_idxlist_delete(gather_list)
    gather_redist = xt_redist_p2p_ext_new(gather_xmap, &
         RESHAPE(scatter_exts, (/ num_parts /)), gather_ext, &
         mpi_double_precision)
    CALL xt_xmap_delete(gather_xmap)
    CALL xt_redist_s_exchange(gather_redist, scattered_data, gathered_data)
    IF (check_correctness_of_exchange) THEN
      DO k = 1, gshape(3)
        DO j = 1, gshape(2)
          DO i = 1, gshape(1)
            IF (gathered_data(i, j, k) < 1.0d0 &
                 .OR. gathered_data(i, j, k) > DBLE(num_parts)) THEN
              WRITE (0, "(a,3(i0,a),g28.16)") "gathered_data(", &
                   i, ", ", j, ", ", k, ") = ", gathered_data(i, j, k)
              WRITE (msg, "(a,3(', ',i0))") "error in data", i, j, k
              CALL xt_abort(msg, filename, __LINE__)
            END IF
          END DO
        END DO
      END DO
    END IF
    CALL xt_redist_delete(gather_redist)
  END DO
  CALL xt_finalize
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) THEN
    WRITE (0, *) "MPI finalization failed"
    CALL posix_exit(1)
  END IF

CONTAINS
  SUBROUTINE parse_options
    INTEGER :: i, num_cmd_args, arg_len
    INTEGER, PARAMETER :: max_opt_arg_len = 80
    CHARACTER(max_opt_arg_len) :: optarg
    CHARACTER(len=18), PARAMETER :: &
         check_correctness_argstr = '-check-correctness'
    num_cmd_args = COMMAND_ARGUMENT_COUNT()
    i = 1
    DO WHILE (i < num_cmd_args)
      CALL GET_COMMAND_ARGUMENT(i, optarg, arg_len)
      IF (optarg(1:arg_len) == check_correctness_argstr &
           .AND. arg_len == LEN(check_correctness_argstr)) THEN
        check_correctness_of_exchange = .TRUE.
      ELSE
        WRITE (0, *) 'unexpected command-line argument parsing error: ', &
             optarg(1:arg_len)
        FLUSH(0)
        CALL xt_abort('unexpected command-line argument "' &
             // optarg(1:arg_len) // '"', filename, __LINE__)
      END IF
      i = i + 1
    END DO

  END SUBROUTINE parse_options

END PROGRAM perf_sparse_array_gather
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
! license-default: "bsd"
! End:
!
