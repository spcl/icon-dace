!>
!! @file xt_idxlist_f.f90
!! @brief Fortran interface to yaxt idxlist methods
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
MODULE xt_idxlist_abstract
  USE xt_core, ONLY: xt_int_kind, xt_mpi_fint_kind, xt_stripe, &
       xt_bounds, xt_abort, i2, i4, i8, &
       xt_pos_ext, OPERATOR(==)
  USE iso_c_binding, ONLY: c_ptr, c_int, c_f_pointer, c_null_ptr, &
       c_associated, c_loc
  IMPLICIT NONE
  PRIVATE

  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_idxlist
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_idxlist

  INTERFACE

    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_idxlist_f2c(idxlist) BIND(c, name='xt_idxlist_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_idxlist
      IMPLICIT NONE
      TYPE(xt_idxlist), INTENT(in) :: idxlist
      TYPE(c_ptr) :: p
    END FUNCTION xt_idxlist_f2c

    FUNCTION xt_idxlist_get_pack_size(idxlist, comm) &
         BIND(c, name='xt_idxlist_get_pack_size_f2c') RESULT(pack_size)
      IMPORT :: xt_idxlist, xt_mpi_fint_kind
      IMPLICIT NONE
      TYPE(xt_idxlist), INTENT(in) :: idxlist
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
      INTEGER(xt_mpi_fint_kind) :: pack_size
    END FUNCTION xt_idxlist_get_pack_size

  END INTERFACE
  ! xt_idxlist_pack_f(type(xt_idxlist), intent(out) :: idxlist,
  !                  type<*>, intent(inout) :: buffer, integer,
  !                  intent(in) :: buffer_size,
  !                  integer, intent(inout) :: position, integer,
  !                  intent(in) :: comm)
  EXTERNAL :: xt_idxlist_pack_f

  ! xt_idxlist_unpack_f(type(xt_idxlist), intent(out) :: idxlist,
  !                     type<*>, intent(in) :: buffer,
  !                     integer, intent(in) :: buffer_size,
  !                     integer, intent(inout) :: position, integer,
  !                     intent(in) :: comm)
  EXTERNAL :: xt_idxlist_unpack_f

  PUBLIC :: xt_idxlist_delete, xt_idxlist_get_pack_size, &
       xt_idxlist_f2c, xt_idxlist_c2f, xt_is_null, &
       xt_idxlist_pack_f, xt_idxlist_unpack_f, xt_idxlist_copy, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_idxlist_get_positions_of_indices, xt_idxlist_get_index_stripes, &
       xt_idxlist_get_num_indices, &
       xt_idxlist_get_indices, xt_idxlist_get_indices_const, &
       xt_idxlist_get_bounding_box, xt_idxlist_get_intersection, &
       xt_idxlist_get_pos_exts_of_index_stripes
  INTERFACE

    FUNCTION xt_idxlist_get_num_indices_c(idxlist) RESULT(num_indices) &
         BIND(c, name='xt_idxlist_get_num_indices')
      IMPORT :: c_int, c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
      INTEGER(c_int) :: num_indices
    END FUNCTION xt_idxlist_get_num_indices_c

    SUBROUTINE xt_idxlist_get_indices_c(idxlist, indices) &
         BIND(c, name='xt_idxlist_get_indices')
      IMPORT :: c_ptr, xt_int_kind
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
      INTEGER(xt_int_kind), INTENT(out) :: indices(*)
    END SUBROUTINE xt_idxlist_get_indices_c

    SUBROUTINE xt_idxlist_delete_c(idxlist) BIND(C, name='xt_idxlist_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
    END SUBROUTINE xt_idxlist_delete_c

    FUNCTION xt_idxlist_get_indices_at_positions_c(idxlist, positions, &
         num_pos, indices, undef_idx) &
         BIND(c, name='xt_idxlist_get_indices_at_positions') RESULT(num_subst)
      IMPORT :: c_ptr, c_int, xt_int_kind
      TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
      INTEGER(c_int), INTENT(in) :: positions(*)
      INTEGER(c_int), VALUE, INTENT(in) :: num_pos
      INTEGER(xt_int_kind), VALUE, INTENT(in) :: undef_idx
      INTEGER(xt_int_kind), INTENT(out) :: indices(*)
      INTEGER(c_int) :: num_subst
    END FUNCTION xt_idxlist_get_indices_at_positions_c

  END INTERFACE

  INTERFACE xt_idxlist_delete
    MODULE PROCEDURE xt_idxlist_delete_1
    MODULE PROCEDURE xt_idxlist_delete_a1d
    MODULE PROCEDURE xt_idxlist_delete_a2d
  END INTERFACE xt_idxlist_delete

  INTERFACE xt_idxlist_get_indices
    MODULE PROCEDURE xt_idxlist_get_indices_1d
    MODULE PROCEDURE xt_idxlist_get_indices_2d
    MODULE PROCEDURE xt_idxlist_get_indices_3d
    MODULE PROCEDURE xt_idxlist_get_indices_4d
    MODULE PROCEDURE xt_idxlist_get_indices_5d
    MODULE PROCEDURE xt_idxlist_get_indices_6d
    MODULE PROCEDURE xt_idxlist_get_indices_7d
  END INTERFACE xt_idxlist_get_indices

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_idxlist_is_null
  END INTERFACE xt_is_null

  INTERFACE xt_idxlist_get_indices_at_positions
    MODULE PROCEDURE xt_idxlist_get_indices_at_positions_a1d
    MODULE PROCEDURE xt_idxlist_get_indices_at_positions_a1d_i2
    MODULE PROCEDURE xt_idxlist_get_indices_at_positions_a1d_i4
    MODULE PROCEDURE xt_idxlist_get_indices_at_positions_a1d_i8
  END INTERFACE xt_idxlist_get_indices_at_positions

  INTERFACE
    FUNCTION xt_idxlist_get_pos_exts_of_index_stripes_c(idxlist, &
       num_stripes, stripes, num_ext, pos_ext, single_match_only) &
       BIND(c, name='xt_idxlist_get_pos_exts_of_index_stripes') &
       RESULT(num_unmatched)
      IMPORT :: c_ptr, c_int
      TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
      INTEGER(c_int), VALUE, INTENT(in) ::  num_stripes
      TYPE(c_ptr), VALUE, INTENT(in) :: stripes
      INTEGER(c_int), INTENT(out) :: num_ext
      TYPE(c_ptr), INTENT(out) :: pos_ext
      INTEGER(c_int), VALUE, INTENT(in) :: single_match_only
      INTEGER(c_int) :: num_unmatched
    END FUNCTION xt_idxlist_get_pos_exts_of_index_stripes_c

    SUBROUTINE free_c(p) BIND(c, name='free')
      IMPORT :: c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: p
    END SUBROUTINE free_c
  END INTERFACE

  INTERFACE xt_idxlist_get_pos_exts_of_index_stripes
    MODULE PROCEDURE gpe_is_i4_a_i4_p1d_l
    MODULE PROCEDURE gpe_is_a_p1d_l
  END INTERFACE xt_idxlist_get_pos_exts_of_index_stripes

  CHARACTER(len=*), PARAMETER :: filename = 'xt_idxlist_f.f90'
CONTAINS

  FUNCTION xt_idxlist_is_null(idxlist) RESULT(p)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(idxlist%cptr)
  END FUNCTION xt_idxlist_is_null

  SUBROUTINE xt_idxlist_delete_1(idxlist)
    TYPE(xt_idxlist), INTENT(inout) :: idxlist
    CALL xt_idxlist_delete_c(idxlist%cptr)
    idxlist%cptr = c_null_ptr
  END SUBROUTINE xt_idxlist_delete_1

  SUBROUTINE xt_idxlist_delete_a1d(idxlists)
    TYPE(xt_idxlist), INTENT(inout) :: idxlists(:)
    INTEGER :: i, n
    n = SIZE(idxlists)
    DO i = 1, n
      CALL xt_idxlist_delete_c(idxlists(i)%cptr)
      idxlists(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_idxlist_delete_a1d

  SUBROUTINE xt_idxlist_delete_a2d(idxlists)
    TYPE(xt_idxlist), INTENT(inout) :: idxlists(:, :)
    INTEGER :: i, j, m, n
    m = SIZE(idxlists, 1)
    n = SIZE(idxlists, 2)
    DO j = 1, n
      DO i = 1, m
        CALL xt_idxlist_delete_c(idxlists(i, j)%cptr)
        idxlists(i, j)%cptr = c_null_ptr
      END DO
    END DO
  END SUBROUTINE xt_idxlist_delete_a2d

  FUNCTION xt_idxlist_get_index_at_position(idxlist, position, idx) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER, VALUE, INTENT(in) :: position
    INTEGER(xt_int_kind), INTENT(out) :: idx
    LOGICAL :: res
    INTEGER(c_int) :: position_c

    INTERFACE
      FUNCTION xt_idxlist_get_index_at_position_c(idxlist, position, idx) &
           BIND(c, name='xt_idxlist_get_index_at_position') RESULT(res)
        IMPORT :: c_ptr, c_int, xt_int_kind
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(c_int), VALUE, INTENT(in) :: position
        INTEGER(xt_int_kind), INTENT(out) :: idx
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_index_at_position_c
    END INTERFACE

    position_c = INT(position, c_int)
    res = xt_idxlist_get_index_at_position_c(idxlist%cptr, &
         position_c, idx) /= 0
  END FUNCTION xt_idxlist_get_index_at_position

  FUNCTION xt_idxlist_get_indices_at_positions_a1d(idxlist, positions, &
    indices, undef_idx) RESULT(num_subst)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER, INTENT(in) :: positions(:)
    INTEGER(xt_int_kind), INTENT(out) :: indices(:)
    INTEGER(xt_int_kind), INTENT(in) :: undef_idx
    INTEGER :: num_subst, n
    INTEGER(c_int) :: num_positions_c

    n = SIZE(positions)
    IF (n > HUGE(1_c_int)) n = HUGE(1_c_int)

    num_positions_c = INT(n, c_int)
    num_subst = xt_idxlist_get_indices_at_positions_c(idxlist%cptr, &
         INT(positions, c_int), num_positions_c, &
         indices, undef_idx)
  END FUNCTION xt_idxlist_get_indices_at_positions_a1d

  FUNCTION xt_idxlist_get_indices_at_positions_a1d_i2(idxlist, positions, &
       num_pos, indices, undef_idx) RESULT(num_subst)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER, INTENT(in) :: positions(*)
    INTEGER(xt_int_kind), INTENT(out) :: indices(*)
    INTEGER(xt_int_kind), INTENT(in) :: undef_idx
    INTEGER(i2), INTENT(in) :: num_pos
    INTEGER :: num_subst
    INTEGER(c_int) :: num_pos_c

    num_pos_c = INT(num_pos, c_int)
    num_subst = xt_idxlist_get_indices_at_positions_c(idxlist%cptr, &
         INT(positions(1:num_pos), c_int), num_pos_c, &
         indices, undef_idx)
  END FUNCTION xt_idxlist_get_indices_at_positions_a1d_i2

  FUNCTION xt_idxlist_get_indices_at_positions_a1d_i4(idxlist, positions, &
       num_pos, indices, undef_idx) RESULT(num_subst)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER, INTENT(in) :: positions(*)
    INTEGER(xt_int_kind), INTENT(out) :: indices(*)
    INTEGER(xt_int_kind), INTENT(in) :: undef_idx
    INTEGER(i4), INTENT(in) :: num_pos
    INTEGER :: num_subst
    INTEGER(c_int) :: num_pos_c

    IF (num_pos > HUGE(1_c_int) .OR. num_pos < 0) &
         CALL xt_abort("invalid number of positions", filename, __LINE__)

    num_pos_c = INT(num_pos, c_int)
    num_subst = xt_idxlist_get_indices_at_positions_c(idxlist%cptr, &
         INT(positions(1:num_pos), c_int), num_pos_c, indices, undef_idx)
  END FUNCTION xt_idxlist_get_indices_at_positions_a1d_i4

  FUNCTION xt_idxlist_get_indices_at_positions_a1d_i8(idxlist, positions, &
       num_pos, indices, undef_idx) RESULT(num_subst)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER, INTENT(in) :: positions(*)
    INTEGER(xt_int_kind), INTENT(out) :: indices(*)
    INTEGER(xt_int_kind), INTENT(in) :: undef_idx
    INTEGER(i8), INTENT(in) :: num_pos
    INTEGER :: num_subst
    INTEGER(c_int) :: num_pos_c

    IF (num_pos > HUGE(1_c_int) .OR. num_pos < 0) &
         CALL xt_abort("invalid number of positions", filename, __LINE__)

    num_pos_c = INT(num_pos, c_int)
    num_subst = xt_idxlist_get_indices_at_positions_c(idxlist%cptr, &
         INT(positions(1:num_pos), c_int), num_pos_c, indices, undef_idx)
  END FUNCTION xt_idxlist_get_indices_at_positions_a1d_i8

  FUNCTION xt_idxlist_get_position_of_index(idxlist, idx, position) &
       RESULT(notfound)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
    INTEGER, INTENT(out) :: position
    LOGICAL :: notfound
    INTEGER(c_int) :: position_c

    INTERFACE
      FUNCTION xt_idxlist_get_position_of_index_c(idxlist, idx, position) &
           BIND(c, name='xt_idxlist_get_position_of_index') RESULT(res)
        IMPORT :: Xt_idxlist, xt_int_kind, c_int, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
        INTEGER(c_int), INTENT(out) :: position
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_position_of_index_c
    END INTERFACE

    notfound = xt_idxlist_get_position_of_index_c(idxlist%cptr, &
         idx, position_c) /= 0
    position = INT(position_c)
  END FUNCTION xt_idxlist_get_position_of_index

  FUNCTION xt_idxlist_get_position_of_index_off(idxlist, idx, position, &
       offset) RESULT(notfound)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
    INTEGER, INTENT(out) :: position
    INTEGER, INTENT(in) :: offset
    LOGICAL :: notfound
    INTEGER(c_int) :: position_c, offset_c

    INTERFACE
      FUNCTION xt_idxlist_get_position_of_index_off_c(idxlist, idx, position, &
           offset) BIND(c, name='xt_idxlist_get_position_of_index_off') &
           RESULT(res)
        IMPORT :: Xt_idxlist, xt_int_kind, c_int, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: idx
        INTEGER(c_int), INTENT(out) :: position
        INTEGER(c_int), VALUE, INTENT(in) :: offset
        INTEGER(c_int) :: res
      END FUNCTION xt_idxlist_get_position_of_index_off_c
    END INTERFACE

    offset_c = INT(offset, c_int)
    notfound = xt_idxlist_get_position_of_index_off_c(idxlist%cptr, &
         idx, position_c, offset_c) /= 0
    position = INT(position_c)
  END FUNCTION xt_idxlist_get_position_of_index_off

  FUNCTION xt_idxlist_get_positions_of_indices(idxlist, indices, positions, &
       single_match_only) RESULT(num_missing)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: indices(:)
    INTEGER, INTENT(out) :: positions(:)
    LOGICAL, INTENT(in) :: single_match_only
    INTEGER :: num_missing, n, ofs
    INTEGER(c_int) :: single_match_only_, num_pos_c

    INTERFACE
      FUNCTION xt_idxlist_get_positions_of_indices_c(idxlist, indices, &
           num_indices, positions, single_match_only) &
           BIND(c, name='xt_idxlist_get_positions_of_indices') &
           RESULT(num_missing)
        IMPORT :: Xt_idxlist, xt_int_kind, c_int, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(xt_int_kind), INTENT(in) :: indices(*)
        INTEGER(c_int), VALUE, INTENT(in) :: num_indices
        INTEGER(c_int), INTENT(out) :: positions(*)
        INTEGER(c_int), VALUE, INTENT(in) :: single_match_only
        INTEGER(c_int) :: num_missing
      END FUNCTION xt_idxlist_get_positions_of_indices_c
    END INTERFACE

    n = SIZE(indices)
    IF (SIZE(positions) < n) THEN
      CALL xt_abort("positions array too small", filename, __LINE__)
    END IF
    num_missing = 0
    ofs = 1
    single_match_only_ = MERGE(1_c_int, 0_c_int, single_match_only)
    DO WHILE (n > 0)
      IF (n > HUGE(1_c_int)) THEN
        num_missing = num_missing &
             + INT(xt_idxlist_get_positions_of_indices_c(&
             idxlist%cptr, indices(ofs:), HUGE(1_c_int), &
             positions(ofs:), single_match_only_))
        ofs = ofs + HUGE(1_c_int)
        n = n - HUGE(1_c_int)
      ELSE
        num_pos_c = INT(n, c_int)
        num_missing = num_missing &
             + INT(xt_idxlist_get_positions_of_indices_c(&
             idxlist%cptr, indices(ofs:), &
             num_pos_c, positions(ofs:), single_match_only_))
        n = 0
      END IF
    END DO
  END FUNCTION xt_idxlist_get_positions_of_indices

  SUBROUTINE xt_idxlist_get_index_stripes(idxlist, stripes)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_stripe), ALLOCATABLE, INTENT(out) :: stripes(:)

    INTERFACE
      SUBROUTINE xt_idxlist_get_index_stripes_c(idxlist, stripes,&
           num_stripes) BIND(c, name='xt_idxlist_get_index_stripes')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        TYPE(c_ptr), INTENT(out) :: stripes
        INTEGER(c_int), INTENT(out) :: num_stripes
      END SUBROUTINE xt_idxlist_get_index_stripes_c
    END INTERFACE
    TYPE(c_ptr) :: stripes_c_ptr
    INTEGER(c_int) :: num_stripes
    TYPE(xt_stripe), POINTER :: stripes_f_ptr(:)
    INTEGER :: stripes_shape(1)
    CALL xt_idxlist_get_index_stripes_c(idxlist%cptr, &
         stripes_c_ptr, num_stripes)
    IF (num_stripes > HUGE(stripes_shape)) &
         CALL xt_abort("number of stripes too large", filename, __LINE__)
    stripes_shape(1) = INT(num_stripes)
    IF (num_stripes > 0) THEN
      ALLOCATE(stripes(INT(num_stripes)))
      CALL C_F_POINTER(stripes_c_ptr, stripes_f_ptr, stripes_shape)
      stripes = stripes_f_ptr
    END IF
    CALL free_c(stripes_c_ptr)
  END SUBROUTINE xt_idxlist_get_index_stripes

  FUNCTION xt_idxlist_get_bounding_box(idxlist, global_size, &
       global_start_index) RESULT(bounds)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(in) :: global_size(:)
    INTEGER(xt_int_kind), INTENT(in) :: global_start_index
    TYPE(xt_bounds) :: bounds(SIZE(global_size))
    INTEGER(c_int) :: ndim

    INTERFACE
      SUBROUTINE xt_idxlist_get_bounding_box_c(idxlist, ndim, global_size, &
           global_start_index, bounds) &
           BIND(c, name='xt_idxlist_get_bounding_box')
        IMPORT :: c_int, c_ptr, xt_int_kind, xt_bounds
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        INTEGER(c_int), VALUE :: ndim
        INTEGER(xt_int_kind), INTENT(in) :: global_size(ndim)
        INTEGER(xt_int_kind), VALUE, INTENT(in) :: global_start_index
        TYPE(xt_bounds), INTENT(out) :: bounds(ndim)
      END SUBROUTINE xt_idxlist_get_bounding_box_c
    END INTERFACE

    ndim = INT(SIZE(global_size), c_int)
    CALL xt_idxlist_get_bounding_box_c(idxlist%cptr, &
         ndim, global_size, global_start_index, bounds)
  END FUNCTION xt_idxlist_get_bounding_box

  FUNCTION xt_idxlist_get_intersection(idxlist_src, idxlist_dst) &
       RESULT(intersection)
    TYPE(xt_idxlist), INTENT(in) :: idxlist_src, idxlist_dst
    TYPE(xt_idxlist) :: intersection

    INTERFACE
      FUNCTION xt_idxlist_get_intersection_c(idxlist_src, idxlist_dst) &
           BIND(c, name='xt_idxlist_get_intersection') RESULT(intersection)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist_src, idxlist_dst
        TYPE(c_ptr) :: intersection
      END FUNCTION xt_idxlist_get_intersection_c
    END INTERFACE

    intersection = xt_idxlist_c2f(xt_idxlist_get_intersection_c(&
         idxlist_src%cptr, idxlist_dst%cptr))
  END FUNCTION xt_idxlist_get_intersection

  FUNCTION xt_idxlist_copy(idxlist) RESULT(copy)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_idxlist) :: copy

    INTERFACE
      FUNCTION xt_idxlist_copy_c(idxlist) BIND(c, name='xt_idxlist_copy') &
           RESULT(copy)
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(IN) :: idxlist
        TYPE(c_ptr) :: copy
      END FUNCTION xt_idxlist_copy_c
    END INTERFACE

    copy = xt_idxlist_c2f(xt_idxlist_copy_c(idxlist%cptr))

  END FUNCTION xt_idxlist_copy

  FUNCTION xt_idxlist_c2f(idxlist) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: idxlist
    TYPE(xt_idxlist) :: p
    p%cptr = idxlist
  END FUNCTION xt_idxlist_c2f

  FUNCTION xt_idxlist_get_num_indices(idxlist) RESULT(num_indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER :: num_indices
    INTEGER(c_int) :: n
    n = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (n > HUGE(num_indices) .OR. n < -HUGE(num_indices)) &
         CALL xt_abort("num_indices out of bounds", filename, __LINE__)
    num_indices = INT(n)
  END FUNCTION xt_idxlist_get_num_indices

  SUBROUTINE xt_idxlist_get_indices_1d(idxlist, indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(out) :: indices(:)
    INTEGER(c_int) :: num_indices
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > SIZE(indices)) THEN
      CALL xt_abort("indices array too small", filename, __LINE__)
    END IF
    CALL xt_idxlist_get_indices_c(idxlist%cptr, indices)
  END SUBROUTINE xt_idxlist_get_indices_1d

  SUBROUTINE xt_idxlist_get_indices_2d(idxlist, indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(out) :: indices(:,:)
    INTEGER(c_int) :: num_indices
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > SIZE(indices)) THEN
      CALL xt_abort("indices array too small", filename, __LINE__)
    END IF
    CALL xt_idxlist_get_indices_c(idxlist%cptr, indices)
  END SUBROUTINE xt_idxlist_get_indices_2d

  SUBROUTINE xt_idxlist_get_indices_3d(idxlist, indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(out) :: indices(:,:,:)
    INTEGER(c_int) :: num_indices
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > SIZE(indices)) THEN
      CALL xt_abort("indices array too small", filename, __LINE__)
    END IF
    CALL xt_idxlist_get_indices_c(idxlist%cptr, indices)
  END SUBROUTINE xt_idxlist_get_indices_3d

  SUBROUTINE xt_idxlist_get_indices_4d(idxlist, indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(out) :: indices(:,:,:,:)
    INTEGER(c_int) :: num_indices
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > SIZE(indices)) THEN
      CALL xt_abort("indices array too small", filename, __LINE__)
    END IF
    CALL xt_idxlist_get_indices_c(idxlist%cptr, indices)
  END SUBROUTINE xt_idxlist_get_indices_4d

  SUBROUTINE xt_idxlist_get_indices_5d(idxlist, indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(out) :: indices(:,:,:,:,:)
    INTEGER(c_int) :: num_indices
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > SIZE(indices)) THEN
      CALL xt_abort("indices array too small", filename, __LINE__)
    END IF
    CALL xt_idxlist_get_indices_c(idxlist%cptr, indices)
  END SUBROUTINE xt_idxlist_get_indices_5d

  SUBROUTINE xt_idxlist_get_indices_6d(idxlist, indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(out) :: indices(:,:,:,:,:,:)
    INTEGER(c_int) :: num_indices
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > SIZE(indices)) THEN
      CALL xt_abort("indices array too small", filename, __LINE__)
    END IF
    CALL xt_idxlist_get_indices_c(idxlist%cptr, indices)
  END SUBROUTINE xt_idxlist_get_indices_6d

  SUBROUTINE xt_idxlist_get_indices_7d(idxlist, indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), INTENT(out) :: indices(:,:,:,:,:,:,:)
    INTEGER(c_int) :: num_indices
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > SIZE(indices)) THEN
      CALL xt_abort("indices array too small", filename, __LINE__)
    END IF
    CALL xt_idxlist_get_indices_c(idxlist%cptr, indices)
  END SUBROUTINE xt_idxlist_get_indices_7d

  FUNCTION xt_idxlist_get_indices_const(idxlist) RESULT(indices)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(xt_int_kind), POINTER :: indices(:)
    INTEGER(c_int) :: num_indices
    TYPE(c_ptr) :: c_indices
    INTEGER(xt_int_kind), SAVE, TARGET :: dummy(1) = -HUGE(indices)
    INTEGER :: indices_shape(1)
    INTERFACE
      FUNCTION xt_idxlist_get_indices_const_c(idxlist) &
           BIND(c, name='xt_idxlist_get_indices_const') RESULT(indices)
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: idxlist
        TYPE(c_ptr) :: indices
      END FUNCTION xt_idxlist_get_indices_const_c
    END INTERFACE
    num_indices = xt_idxlist_get_num_indices_c(idxlist%cptr)
    IF (num_indices > 0_xt_int_kind) THEN
      IF (num_indices > HUGE(indices_shape)) &
           CALL xt_abort("too many indices for default integer kind", &
           filename, __LINE__)
      indices_shape(1) = INT(num_indices)
      c_indices = xt_idxlist_get_indices_const_c(idxlist%cptr)
      CALL C_F_POINTER(c_indices, indices, indices_shape)
    ELSE
      indices => dummy(1:0)
    END IF
  END FUNCTION xt_idxlist_get_indices_const

  FUNCTION gpe_is_i4_a_i4_p1d_l(idxlist, &
       num_stripes, stripes, num_ext, pos_ext, single_match_only) &
       RESULT(num_unmatched)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    INTEGER(i4), INTENT(in) ::  num_stripes
    TYPE(xt_stripe), INTENT(in), TARGET :: stripes(num_stripes)
    INTEGER, INTENT(out) :: num_ext
    TYPE(xt_pos_ext), ALLOCATABLE, INTENT(out) :: pos_ext(:)
    LOGICAL, INTENT(in) :: single_match_only
    INTEGER :: num_unmatched

    INTEGER(c_int) :: num_unmatched_c, num_ext_c, num_stripes_c
    TYPE(c_ptr) :: pos_ext_c, stripes_c
    TYPE(xt_pos_ext), POINTER :: pos_ext_fptr(:)
    INTEGER :: pos_ext_shape(1)
    TYPE(xt_pos_ext), TARGET :: dummy_stripe(1)

    IF (num_stripes > HUGE(1_c_int) .OR. num_stripes < 0) &
         CALL xt_abort("interface violation detected", filename, __LINE__)

    IF (num_stripes > 0_i4) THEN
      stripes_c = C_LOC(stripes)
    ELSE
      stripes_c = C_LOC(dummy_stripe)
    END IF
    num_stripes_c = INT(num_stripes, c_int)
    num_unmatched_c = xt_idxlist_get_pos_exts_of_index_stripes_c(&
         idxlist%cptr, num_stripes_c, stripes_c, &
         num_ext_c, pos_ext_c, MERGE(1_c_int, 0_c_int, single_match_only))

    IF (num_ext_c > HUGE(1) .OR. num_ext_c < 0 &
         .OR. num_unmatched_c > HUGE(1) .OR. num_unmatched_c < 0) &
         CALL xt_abort("data representation problem", filename, __LINE__)
    num_unmatched = INT(num_unmatched_c)
    num_ext = INT(num_ext_c)
    IF (num_ext > 0) THEN
      ALLOCATE(pos_ext(num_ext))
      pos_ext_shape(1) = num_ext
      CALL C_F_POINTER(pos_ext_c, pos_ext_fptr, pos_ext_shape)
      pos_ext = pos_ext_fptr
      CALL free_c(pos_ext_c)
    END IF
  END FUNCTION gpe_is_i4_a_i4_p1d_l

  FUNCTION gpe_is_a_p1d_l(idxlist, stripes, pos_ext, single_match_only) &
       RESULT(num_unmatched)
    TYPE(xt_idxlist), INTENT(in) :: idxlist
    TYPE(xt_stripe), INTENT(in) :: stripes(:)
    TYPE(xt_pos_ext), ALLOCATABLE, INTENT(out) :: pos_ext(:)
    LOGICAL, INTENT(in) :: single_match_only
    INTEGER :: num_unmatched

    INTEGER :: num_ext
    INTEGER(i4) :: num_stripes

    num_stripes = SIZE(stripes)
    IF (num_stripes > 0) THEN
      num_unmatched = gpe_is_i4_a_i4_p1d_l(idxlist, num_stripes, stripes, &
           num_ext, pos_ext, single_match_only)
    ELSE
      num_unmatched = 0
    END IF
  END FUNCTION gpe_is_a_p1d_l

END MODULE xt_idxlist_abstract
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
