!>
!! @file xt_sort_f.f90
!! @brief Fortran interface to yaxt sort declarations
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
MODULE xt_sort
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_size_t
  USE xt_core, ONLY: xt_abort, xt_int_kind
  IMPLICIT NONE
  PRIVATE

  INTERFACE xt_sort_int
    SUBROUTINE xt_sort_int_f2c(a, n) BIND(c, name='xt_sort_int_f2c')
      IMPORT :: c_int, c_size_t
      INTEGER(c_size_t), VALUE, INTENT(in) :: n
      INTEGER(c_int), INTENT(inout) :: a(n)
    END SUBROUTINE xt_sort_int_f2c
    MODULE PROCEDURE xt_sort_int_a
  END INTERFACE xt_sort_int
  PUBLIC :: xt_sort_int

  INTERFACE xt_sort_index
    SUBROUTINE xt_sort_index_f2c(a, n, positions, reset_positions) &
          BIND(c, name='xt_sort_index_f2c')
      IMPORT :: c_int, xt_int_kind
      INTEGER(c_int), VALUE, INTENT(in) :: n, reset_positions
      INTEGER(xt_int_kind), INTENT(inout) :: a(n)
      INTEGER(c_int), INTENT(inout) :: positions(n)
    END SUBROUTINE xt_sort_index_f2c
    MODULE PROCEDURE xt_sort_index_a_a_l
  END INTERFACE xt_sort_index
  PUBLIC :: xt_sort_index

  INTERFACE xt_sort_permutation
    SUBROUTINE xt_sort_int_permutation_f2c(a, n, permutation) &
          BIND(c, name='xt_sort_int_permutation_f2c')
      IMPORT :: c_int, c_size_t
      INTEGER(c_size_t), VALUE, INTENT(in) :: n
      INTEGER(c_int), INTENT(inout) :: a(n), permutation(n)
    END SUBROUTINE xt_sort_int_permutation_f2c
    MODULE PROCEDURE xt_sort_int_permutation
  END INTERFACE xt_sort_permutation
  PUBLIC :: xt_sort_permutation

  TYPE, BIND(c), PUBLIC :: xt_idxpos
    INTEGER(xt_int_kind) :: idx
    INTEGER(c_int) :: pos
  END TYPE xt_idxpos

  INTERFACE xt_sort_idxpos
    SUBROUTINE xt_sort_idxpos_f2c(a, n) BIND(c, name='xt_sort_idxpos_f2c')
      IMPORT :: c_size_t, xt_idxpos
      INTEGER(c_size_t), VALUE, INTENT(in) :: n
      TYPE(xt_idxpos), INTENT(inout) :: a(n)
    END SUBROUTINE xt_sort_idxpos_f2c
    MODULE PROCEDURE xt_sort_idxpos_a
  END INTERFACE xt_sort_idxpos
  PUBLIC :: xt_sort_idxpos

  INTERFACE xt_assign_id_map
    SUBROUTINE xt_assign_id_map_int(n, a, ofs) &
         BIND(c, name='xt_assign_id_map_int')
      IMPORT :: c_int, c_size_t
      INTEGER(c_size_t), VALUE, INTENT(in) :: n
      INTEGER(c_int), INTENT(out) :: a(n)
      INTEGER(c_int), VALUE, INTENT(in) :: ofs
    END SUBROUTINE xt_assign_id_map_int
    MODULE PROCEDURE xt_assign_id_map_a
    MODULE PROCEDURE xt_assign_id_map_ai
  END INTERFACE xt_assign_id_map
  PUBLIC :: xt_assign_id_map

  CHARACTER(len=*), PARAMETER :: filename = 'xt_sort_f.f90'
CONTAINS
  SUBROUTINE xt_sort_int_a(a)
    INTEGER(c_int), INTENT(inout) :: a(:)
    INTEGER :: a_size
    a_size = SIZE(a)
    IF (a_size > 1) CALL xt_sort_int_f2c(a, INT(a_size, c_size_t))
  END SUBROUTINE xt_sort_int_a

  SUBROUTINE xt_sort_index_a_a_l(a, positions, reset_positions)
    INTEGER(xt_int_kind), INTENT(inout) :: a(:)
    INTEGER(c_int), INTENT(inout) :: positions(:)
    LOGICAL, INTENT(in) :: reset_positions
    INTEGER(c_int) :: a_size, reset_positions_c
    a_size = INT(SIZE(a), c_int)
    reset_positions_c = MERGE(1_c_int, 0_c_int, reset_positions)
    IF (a_size > SIZE(positions)) &
         CALL xt_abort("positions array too small", filename, __LINE__)
    IF (a_size > 1) CALL xt_sort_index_f2c(a, a_size, positions, &
         reset_positions_c)
  END SUBROUTINE xt_sort_index_a_a_l

  SUBROUTINE xt_sort_idxpos_a(a)
    TYPE(xt_idxpos), INTENT(inout) :: a(:)
    INTEGER :: a_size
    a_size = SIZE(a)
    IF (a_size > 1) CALL xt_sort_idxpos_f2c(a, INT(a_size, c_size_t))
  END SUBROUTINE xt_sort_idxpos_a

  SUBROUTINE xt_sort_int_permutation(a, permutation)
    INTEGER(c_int), INTENT(inout) :: a(:), permutation(:)
    INTEGER(c_size_t) :: a_size
#ifdef HAVE_SIZE_KIND_ARGUMENT
    a_size = SIZE(a, kind=c_size_t)
#else
    a_size = INT(SIZE(a), kind=c_size_t)
#endif
    IF (a_size > 1) CALL xt_sort_int_permutation_f2c(a, a_size, permutation)
  END SUBROUTINE xt_sort_int_permutation

  SUBROUTINE xt_assign_id_map_a(a)
    INTEGER(c_int), INTENT(out) :: a(:)
    INTEGER(c_size_t) :: n
    INTEGER(c_int) :: ofs
#ifdef HAVE_SIZE_KIND_ARGUMENT
    n = SIZE(a, kind=c_size_t)
#else
    n = INT(SIZE(a), kind=c_size_t)
#endif
    IF (n >= 1) THEN
      ofs = 0_c_int
      CALL xt_assign_id_map_int(n, a, ofs)
    END IF
  END SUBROUTINE xt_assign_id_map_a

  SUBROUTINE xt_assign_id_map_ai(a, ofs)
    INTEGER(c_int), INTENT(out) :: a(:)
    INTEGER(c_int), INTENT(in) :: ofs
    INTEGER(c_size_t) :: n
#ifdef HAVE_SIZE_KIND_ARGUMENT
    n = SIZE(a, kind=c_size_t)
#else
    n = INT(SIZE(a), kind=c_size_t)
#endif
    IF (n >= 1) THEN
      CALL xt_assign_id_map_int(n, a, ofs)
    END IF
  END SUBROUTINE xt_assign_id_map_ai

END MODULE xt_sort
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
