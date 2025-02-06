!>
!! @file xt_xmap_intersection_f.f90
!! @brief Fortran interface to yaxt xmap declarations
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

!>
!! @example test_xmap_intersection_parallel_f.f90

#include "fc_feature_defs.inc"
MODULE xt_xmap_intersection
  USE iso_c_binding, ONLY: c_int, c_loc, c_null_ptr, c_ptr
  USE xt_core, ONLY: xt_abort, xt_mpi_fint_kind
  USE xt_idxlist_abstract, ONLY: xt_idxlist, xt_idxlist_f2c
  USE xt_xmap_abstract, ONLY: xt_xmap, xt_xmap_c2f
#include "xt_slice_c_loc.inc"
  IMPLICIT NONE
  PRIVATE

  TYPE, BIND(c), PUBLIC :: xt_com_list
    TYPE(xt_idxlist) :: list
    INTEGER(c_int) :: rank
  END TYPE xt_com_list

  TYPE, PUBLIC :: xt_com_pos
    INTEGER, POINTER :: transfer_pos(:)
    INTEGER :: rank
  END TYPE xt_com_pos

  TYPE, BIND(c) :: xt_com_pos_c
    TYPE(c_ptr) :: transfer_pos
    INTEGER(c_int) :: num_transfer_pos
    INTEGER(c_int) :: rank
  END TYPE xt_com_pos_c

  PUBLIC :: xt_xmap_intersection_new, xt_xmap_intersection_ext_new, &
            xt_xmap_intersection_pos_new

  INTERFACE
    FUNCTION xmi_new_f2c(num_src_intersections, src_com, &
         num_dst_intersections, dst_com, &
         src_idxlist, dst_idxlist, comm) RESULT(xmap) &
         BIND(c, name='xt_xmap_intersection_new_f2c')
      IMPORT :: c_int, c_ptr, xt_mpi_fint_kind
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
           num_dst_intersections
      TYPE(c_ptr), VALUE, INTENT(in) :: src_com, dst_com, src_idxlist, &
           dst_idxlist
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
      TYPE(c_ptr) :: xmap
    END FUNCTION xmi_new_f2c

    FUNCTION xmi_ext_new_f2c(num_src_intersections, src_com, &
         num_dst_intersections, dst_com, &
         src_idxlist, dst_idxlist, comm) RESULT(xmap) &
         BIND(c, name='xt_xmap_intersection_ext_new_f2c')
      IMPORT :: c_int, c_ptr, xt_mpi_fint_kind
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
           num_dst_intersections
      TYPE(c_ptr), VALUE, INTENT(in) :: src_com, dst_com, src_idxlist, &
           dst_idxlist
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
      TYPE(c_ptr) :: xmap
    END FUNCTION xmi_ext_new_f2c

    FUNCTION xmi_pos_new_f2c(num_src_msg, src_com, num_dst_msg, dst_com, comm) &
         RESULT(xmap) &
         BIND(c, name='xt_xmap_intersection_pos_new_f2c')
      IMPORT :: c_int, c_ptr, xt_mpi_fint_kind, xt_com_pos_c
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_msg, num_dst_msg
      TYPE(xt_com_pos_c), INTENT(in) :: src_com(*), dst_com(*)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
      TYPE(c_ptr) :: xmap
    END FUNCTION xmi_pos_new_f2c
  END INTERFACE

  INTERFACE xt_xmap_intersection_new
    MODULE PROCEDURE xmi_new_i_a_i_a
    MODULE PROCEDURE xmi_new_a_a
  END INTERFACE xt_xmap_intersection_new

  INTERFACE xt_xmap_intersection_ext_new
    MODULE PROCEDURE xmi_ext_new_i_a_i_a
    MODULE PROCEDURE xmi_ext_new_a_a
  END INTERFACE xt_xmap_intersection_ext_new

  INTERFACE xt_xmap_intersection_pos_new
    MODULE PROCEDURE xmi_pos_new_a_a
    MODULE PROCEDURE xmi_pos_new_i_a_i_a
  END INTERFACE xt_xmap_intersection_pos_new

  CHARACTER(len=*), PARAMETER :: filename = 'xt_xmap_intersection_f.f90'
CONTAINS
  FUNCTION xmi_new_i_a_i_a(num_src_intersections, src_com, &
       num_dst_intersections, dst_com, &
       src_idxlist, dst_idxlist, comm) RESULT(xmap)
    INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
         num_dst_intersections
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(num_src_intersections), &
         dst_com(num_dst_intersections)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap
    TYPE(c_ptr) :: src_com_p, dst_com_p

    IF (num_src_intersections > 0) THEN
      src_com_p = C_LOC(src_com)
    ELSE
      src_com_p = c_null_ptr
    END IF
    IF (num_dst_intersections > 0) THEN
      dst_com_p = C_LOC(dst_com)
    ELSE
      dst_com_p = c_null_ptr
    END IF
    xmap = xt_xmap_c2f(xmi_new_f2c(&
         num_src_intersections, src_com_p, &
         num_dst_intersections, dst_com_p, &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_new_i_a_i_a

  FUNCTION xmi_new_a_a(src_com, dst_com, src_idxlist, dst_idxlist, comm) &
       RESULT(xmap)
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(:), dst_com(:)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    TYPE(xt_com_list), ALLOCATABLE, TARGET :: src_com_a(:), dst_com_a(:)
    TYPE(c_ptr) :: src_com_p, dst_com_p
    INTEGER(c_int) :: num_src_intersections_c, num_dst_intersections_c
    num_src_intersections_c = INT(SIZE(src_com), c_int)
    num_dst_intersections_c = INT(SIZE(dst_com), c_int)
    CALL com_p_arg(src_com, src_com_a, src_com_p)
    CALL com_p_arg(dst_com, dst_com_a, dst_com_p)

    xmap = xt_xmap_c2f(xmi_new_f2c(num_src_intersections_c, src_com_p, &
         num_dst_intersections_c, dst_com_p, &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_new_a_a

  FUNCTION xmi_ext_new_i_a_i_a(num_src_intersections, src_com, &
       num_dst_intersections, dst_com, &
       src_idxlist, dst_idxlist, comm) RESULT(xmap)
    INTEGER(c_int), VALUE, INTENT(in) :: num_src_intersections, &
         num_dst_intersections
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(num_src_intersections), &
         dst_com(num_dst_intersections)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap
    INTEGER(c_int) :: num_src_intersections_c, num_dst_intersections_c
    num_src_intersections_c = INT(num_src_intersections, c_int)
    num_dst_intersections_c = INT(num_dst_intersections, c_int)

    xmap = xt_xmap_c2f(xmi_ext_new_f2c(&
         num_src_intersections_c, C_LOC(src_com), &
         num_dst_intersections_c, C_LOC(dst_com), &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_ext_new_i_a_i_a

  FUNCTION xmi_ext_new_a_a(src_com, dst_com, src_idxlist, dst_idxlist, comm) &
       RESULT(xmap)
    TYPE(xt_com_list), TARGET, INTENT(in) :: src_com(:), dst_com(:)
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    TYPE(xt_com_list), ALLOCATABLE, TARGET :: src_com_a(:), dst_com_a(:)
    TYPE(c_ptr) :: src_com_p, dst_com_p
    INTEGER(c_int) :: num_src_intersections_c, num_dst_intersections_c
    num_src_intersections_c = INT(SIZE(src_com), c_int)
    num_dst_intersections_c = INT(SIZE(dst_com), c_int)

    CALL com_p_arg(src_com, src_com_a, src_com_p)
    CALL com_p_arg(dst_com, dst_com_a, dst_com_p)

    xmap = xt_xmap_c2f(xmi_ext_new_f2c(&
         num_src_intersections_c, src_com_p, &
         num_dst_intersections_c, dst_com_p, &
         xt_idxlist_f2c(src_idxlist), xt_idxlist_f2c(dst_idxlist), comm))
  END FUNCTION xmi_ext_new_a_a

  FUNCTION xmi_pos_new_i_a_i_a( &
      num_src_msg, src_com, num_dst_msg, dst_com, comm) RESULT(xmap)
    TYPE(xt_com_pos), TARGET, INTENT(in) :: src_com(:), dst_com(:)
    INTEGER, INTENT(in) :: num_src_msg, num_dst_msg
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    INTEGER(c_int) :: num_src_msg_c, num_dst_msg_c
    TYPE(xt_com_pos_c), ALLOCATABLE :: src_com_c(:), dst_com_c(:)
    INTEGER(c_int), TARGET, ALLOCATABLE :: pos_buffer(:)
    INTEGER :: pos_buffer_offset, size_pos_buf

    size_pos_buf = num_pos_copy(num_src_msg, src_com) + &
                    num_pos_copy(num_dst_msg, dst_com)
    ALLOCATE(pos_buffer(size_pos_buf))

    num_src_msg_c = INT(num_src_msg, c_int)
    num_dst_msg_c = INT(num_dst_msg, c_int)

    pos_buffer_offset = 0
    CALL generate_xt_com_pos_c(num_src_msg, src_com, src_com_c, &
         size_pos_buf, pos_buffer, pos_buffer_offset)
    CALL generate_xt_com_pos_c(num_dst_msg, dst_com, dst_com_c, &
         size_pos_buf, pos_buffer, pos_buffer_offset)

    xmap = &
      xt_xmap_c2f(xmi_pos_new_f2c(&
        num_src_msg_c, src_com_c, num_dst_msg_c, dst_com_c, comm))

  END FUNCTION xmi_pos_new_i_a_i_a

  PURE FUNCTION num_pos_copy(num_msg, com_pos) RESULT(total_num_pos)
    INTEGER, INTENT(in) :: num_msg
    TYPE(xt_com_pos), INTENT(in) :: com_pos(:)
    INTEGER :: i
    INTEGER :: total_num_pos
#if defined __PGI && __PGIC__  > 15 && __PGIC__ <= 20
    INTEGER, POINTER :: pos(:)
#endif
    total_num_pos = 0
#ifdef HAVE_FC_IS_CONTIGUOUS
    IF (KIND(com_pos(i)%transfer_pos) == c_int) THEN
      DO i = 1, num_msg
#if defined __PGI && __PGIC__  > 15 && __PGIC__ <= 20
        pos => com_pos(i)%transfer_pos
        IF (.NOT. IS_CONTIGUOUS(pos)) THEN
#else
        IF (.NOT. IS_CONTIGUOUS(com_pos(i)%transfer_pos)) THEN
#endif
          total_num_pos = total_num_pos &
               + SIZE(com_pos(i)%transfer_pos)
        END IF
      END DO
    ELSE
#endif
      DO i = 1, num_msg
        total_num_pos = total_num_pos + SIZE(com_pos(i)%transfer_pos)
      END DO
#ifdef HAVE_FC_IS_CONTIGUOUS
    ENDIF
#endif
  END FUNCTION num_pos_copy

  SUBROUTINE generate_xt_com_pos_c(num_msg, com_pos, com_pos_c, &
       size_pos_buf, pos_buffer, &
       pos_buffer_offset)
    INTEGER, INTENT(in) :: num_msg
    TYPE(xt_com_pos), TARGET, INTENT(in) :: com_pos(:)
    TYPE(xt_com_pos_c), ALLOCATABLE, INTENT(out) :: com_pos_c(:)
    INTEGER, INTENT(in) :: size_pos_buf
    INTEGER(c_int), TARGET, INTENT(inout) :: pos_buffer(size_pos_buf)
    INTEGER, INTENT(inout) :: pos_buffer_offset
    INTEGER :: i, j, num_pos
#if defined __PGI && __PGIC__  > 15 && __PGIC__ <= 20
    INTEGER, POINTER :: pos(:)
#endif
    ALLOCATE(com_pos_c(num_msg))

    DO i = 1, num_msg
      num_pos = SIZE(com_pos(i)%transfer_pos)
#ifdef HAVE_FC_IS_CONTIGUOUS
#  if defined __PGI && __PGIC__ > 15 && __PGIC__ <= 20
      pos => com_pos(i)%transfer_pos
      IF (KIND(1) == c_int .AND. IS_CONTIGUOUS(pos)) THEN
#  else
      IF (KIND(1) == c_int .AND. IS_CONTIGUOUS(com_pos(i)%transfer_pos)) THEN
#  endif
        com_pos_c(i)%transfer_pos = C_LOC(com_pos(i)%transfer_pos(1))
      ELSE
#endif
        DO j = 1, num_pos
          pos_buffer(pos_buffer_offset + j) = &
                 INT(com_pos(i)%transfer_pos(j), c_int)
        END DO
        com_pos_c(i)%transfer_pos = C_LOC(pos_buffer(pos_buffer_offset+1))
        pos_buffer_offset = pos_buffer_offset + num_pos
#ifdef HAVE_FC_IS_CONTIGUOUS
      END IF
#endif
      com_pos_c(i)%num_transfer_pos = INT(num_pos, c_int)
      com_pos_c(i)%rank = INT(com_pos(i)%rank, c_int)
    END DO

  END SUBROUTINE generate_xt_com_pos_c

  FUNCTION xmi_pos_new_a_a(src_com, dst_com, comm) RESULT(xmap)
    TYPE(xt_com_pos), INTENT(in) :: src_com(:), dst_com(:)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    xmap = &
      xmi_pos_new_i_a_i_a(SIZE(src_com), src_com, SIZE(dst_com), dst_com, comm)
  END FUNCTION xmi_pos_new_a_a

  SUBROUTINE com_p_arg(com, com_a, com_p)
    TYPE(xt_com_list), TARGET, INTENT(in) :: com(:)
    TYPE(xt_com_list), TARGET, ALLOCATABLE, INTENT(inout) :: com_a(:)
    TYPE(c_ptr), INTENT(out) :: com_p

    INTEGER :: com_size
    LOGICAL :: com_is_contiguous
#ifndef HAVE_FC_IS_CONTIGUOUS
    INTERFACE
      FUNCTION xt_com_list_contiguous(com_a, com_b) RESULT(p) &
            BIND(c, name='xt_com_list_contiguous')
        IMPORT :: c_int, xt_com_list
        TYPE(xt_com_list), INTENT(in) :: com_a, com_b
        INTEGER(c_int) :: p
      END FUNCTION xt_com_list_contiguous
    END INTERFACE
#endif

    com_size = SIZE(com)
    IF (com_size > HUGE(1_c_int)) &
         CALL xt_abort('invalid size', filename, __LINE__)
    IF (com_size > 0) THEN
      IF (com_size > 1) THEN
#ifdef HAVE_FC_IS_CONTIGUOUS
        com_is_contiguous = IS_CONTIGUOUS(com)
#else
        com_is_contiguous = xt_com_list_contiguous(com(1), com(2)) /= 0
#endif
        IF (com_is_contiguous) THEN
          XT_SLICE_C_LOC(com(1), com_p)
        ELSE
          ALLOCATE(com_a(com_size))
          com_a = com
          com_p = C_LOC(com_a)
        END IF
      ELSE
        XT_SLICE_C_LOC(com(1), com_p)
      END IF
    ELSE
      com_p = c_null_ptr
    END IF
  END SUBROUTINE com_p_arg

END MODULE xt_xmap_intersection
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
