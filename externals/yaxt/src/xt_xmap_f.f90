!>
!! @file xt_xmap_f.f90
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
#include "fc_feature_defs.inc"
MODULE xt_xmap_abstract
  USE iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, &
       c_associated, c_f_pointer, c_loc
  USE xt_core, ONLY: xt_abort, xt_mpi_fint_kind, xt_pos_ext, i2, i4, i8
  USE xt_config_f, ONLY: xt_config
  USE xt_idxlist_abstract, ONLY: xt_idxlist
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: xt_xmap_c2f, xt_xmap_f2c, xt_is_null
  PUBLIC :: xt_xmap_copy, xt_xmap_delete, xt_xmap_get_num_destinations, &
       xt_xmap_get_num_sources, xt_xmap_get_destination_ranks, &
       xt_xmap_get_max_src_pos, xt_xmap_get_max_dst_pos, &
       xt_xmap_get_source_ranks, xt_xmap_reorder, xt_reorder_type_kind, &
       xt_reorder_none, xt_reorder_send_up, xt_reorder_recv_up, &
       xt_xmap_update_positions, xt_xmap_spread
  PUBLIC :: xt_xmap_all2all_new, xt_xmap_all2all_custom_new
  PUBLIC :: xt_xmap_dist_dir_new,  xt_xmap_dist_dir_custom_new
  PUBLIC ::  xt_xmap_dist_dir_intercomm_new,  &
       xt_xmap_dist_dir_intercomm_custom_new
  PUBLIC :: xt_xmap_get_in_iterator, xt_xmap_get_out_iterator, &
       xt_xmap_iterator_next, xt_xmap_iterator_get_rank, &
       xt_xmap_iterator_get_transfer_pos, &
       xt_xmap_iterator_get_transfer_pos_ext, &
       xt_xmap_iterator_get_num_transfer_pos, &
       xt_xmap_iterator_get_num_transfer_pos_ext, &
       xt_xmap_iterator_delete


  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_xmap
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmap

  TYPE, BIND(c), PUBLIC :: xt_xmap_iter
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmap_iter

  ENUM, BIND( C )
    ENUMERATOR :: xt_reorder_none, xt_reorder_send_up, xt_reorder_recv_up
  END ENUM
  INTEGER, PARAMETER :: xt_reorder_type_kind = KIND(xt_reorder_none)

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_xmap_f2c(xmap) BIND(c, name='xt_xmap_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_xmap
      IMPLICIT NONE
      TYPE(xt_xmap), INTENT(in) :: xmap
      TYPE(c_ptr) :: p
    END FUNCTION xt_xmap_f2c

    SUBROUTINE xt_xmap_delete_c(xmap) BIND(C, name='xt_xmap_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: xmap
    END SUBROUTINE xt_xmap_delete_c

  END INTERFACE

  INTERFACE xt_xmap_delete
    MODULE PROCEDURE xt_xmap_delete_1
    MODULE PROCEDURE xt_xmap_delete_a1d
  END INTERFACE xt_xmap_delete

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_xmap_is_null
    MODULE PROCEDURE xt_xmap_iterator_is_null
  END INTERFACE xt_is_null

  INTERFACE
    FUNCTION xt_xmap_iterator_get_num_transfer_pos_c(iter) RESULT(num) &
         BIND(c, name='xt_xmap_iterator_get_num_transfer_pos')
      IMPORT :: c_int, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iter
      INTEGER(c_int) :: num
    END FUNCTION xt_xmap_iterator_get_num_transfer_pos_c

    FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext_c(iter) RESULT(num) &
         BIND(c, name='xt_xmap_iterator_get_num_transfer_pos_ext')
      IMPORT :: c_int, c_ptr
      TYPE(c_ptr), VALUE, INTENT(in) :: iter
      INTEGER(c_int) :: num
    END FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext_c
  END INTERFACE

  INTERFACE xt_xmap_spread
    MODULE PROCEDURE xt_xmap_spread_a1d
    MODULE PROCEDURE xt_xmap_spread_i2_a1d
    MODULE PROCEDURE xt_xmap_spread_i4_a1d
    MODULE PROCEDURE xt_xmap_spread_i8_a1d
  END INTERFACE xt_xmap_spread

  CHARACTER(len=*), PARAMETER :: filename = 'xt_xmap_f.f90'
CONTAINS

  FUNCTION xt_xmap_is_null(xmap) RESULT(p)
    TYPE(xt_xmap), INTENT(in) :: xmap
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(xmap%cptr)
  END FUNCTION xt_xmap_is_null


  FUNCTION xt_xmap_c2f(xmap) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: xmap
    TYPE(xt_xmap) :: p
    p%cptr = xmap
  END FUNCTION xt_xmap_c2f

  FUNCTION xt_xmap_copy(xmap) RESULT(xmap_copy)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_xmap) :: xmap_copy
    INTERFACE
      FUNCTION xt_xmap_copy_c(xmap) BIND(C, name='xt_xmap_copy') RESULT(res_ptr)
        IMPORT :: xt_xmap, c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_xmap_copy_c
    END INTERFACE
    xmap_copy%cptr = xt_xmap_copy_c(xmap%cptr)
  END FUNCTION xt_xmap_copy
  SUBROUTINE xt_xmap_delete_1(xmap)
    TYPE(xt_xmap), INTENT(inout) :: xmap
    CALL xt_xmap_delete_c(xmap%cptr)
    xmap%cptr = c_null_ptr
  END SUBROUTINE xt_xmap_delete_1

  SUBROUTINE xt_xmap_delete_a1d(xmaps)
    TYPE(xt_xmap), INTENT(inout) :: xmaps(:)
    INTEGER :: i, n
    n = SIZE(xmaps)
    DO i = 1, n
      CALL xt_xmap_delete_c(xmaps(i)%cptr)
      xmaps(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_xmap_delete_a1d

  FUNCTION xt_xmap_all2all_new(src_idxlist, dst_idxlist, comm) RESULT(xmap)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    INTERFACE
      FUNCTION xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm) &
           BIND(C, name='xt_xmap_all2all_new_f') RESULT(xmap_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(c_ptr) :: xmap_ptr
      END FUNCTION xt_xmap_all2all_new_f
    END INTERFACE

    xmap%cptr = xt_xmap_all2all_new_f(src_idxlist, dst_idxlist, comm)
  END FUNCTION xt_xmap_all2all_new

  FUNCTION xt_xmap_all2all_custom_new(src_idxlist, dst_idxlist, comm, config) &
       RESULT(xmap)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_xmap) :: xmap

    INTERFACE
      FUNCTION xt_xmap_all2all_custom_new_f(src_idxlist, dst_idxlist, comm, &
           config) BIND(C, name='xt_xmap_all2all_custom_new_f') RESULT(xmap_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, xt_config, c_ptr
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(xt_config), INTENT(in) :: config
        TYPE(c_ptr) :: xmap_ptr
      END FUNCTION xt_xmap_all2all_custom_new_f
    END INTERFACE

    xmap%cptr = xt_xmap_all2all_custom_new_f(src_idxlist, dst_idxlist, &
         comm, config)
  END FUNCTION xt_xmap_all2all_custom_new

  FUNCTION xt_xmap_dist_dir_new(src_idxlist, dst_idxlist, comm) RESULT(xmap)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap

    INTERFACE
      FUNCTION xt_xmap_dist_dir_new_f(src_idxlist, dst_idxlist, comm) &
           BIND(C, name='xt_xmap_dist_dir_new_f') RESULT(xmap_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(c_ptr) :: xmap_ptr
      END FUNCTION xt_xmap_dist_dir_new_f
    END INTERFACE

    xmap%cptr = xt_xmap_dist_dir_new_f(src_idxlist, dst_idxlist, comm)
  END FUNCTION xt_xmap_dist_dir_new

  FUNCTION xt_xmap_dist_dir_custom_new(src_idxlist, dst_idxlist, comm, config) &
       RESULT(xmap)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    TYPE(xt_config), INTENT(in) :: config
    INTEGER, INTENT(in) :: comm
    TYPE(xt_xmap) :: xmap
    INTERFACE
      FUNCTION xt_xmap_dist_dir_custom_new_f(src_idxlist, dst_idxlist, comm, &
           config) BIND(C, name='xt_xmap_dist_dir_custom_new_f') &
           RESULT(xmap_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, xt_config, c_ptr
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm
        TYPE(xt_config), INTENT(in) :: config
        TYPE(c_ptr) :: xmap_ptr
      END FUNCTION xt_xmap_dist_dir_custom_new_f
    END INTERFACE
    xmap%cptr = xt_xmap_dist_dir_custom_new_f(src_idxlist, dst_idxlist, &
         comm, config)
  END FUNCTION xt_xmap_dist_dir_custom_new

  FUNCTION xt_xmap_dist_dir_intercomm_new(src_idxlist, dst_idxlist, &
       inter_comm, intra_comm) RESULT(xmap)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, INTENT(in) :: inter_comm, intra_comm
    TYPE(xt_xmap) :: xmap

    INTERFACE
      FUNCTION xt_xmap_dist_dir_intercomm_new_f(src_idxlist, dst_idxlist, &
           inter_comm, intra_comm) &
           BIND(C, name='xt_xmap_dist_dir_intercomm_new_f') RESULT(xmap_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, c_ptr
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: inter_comm, intra_comm
        TYPE(c_ptr) :: xmap_ptr
      END FUNCTION xt_xmap_dist_dir_intercomm_new_f
    END INTERFACE

    xmap%cptr = xt_xmap_dist_dir_intercomm_new_f(src_idxlist, &
         dst_idxlist, inter_comm, intra_comm)
  END FUNCTION xt_xmap_dist_dir_intercomm_new

  FUNCTION xt_xmap_dist_dir_intercomm_custom_new(src_idxlist, dst_idxlist, &
       inter_comm, intra_comm, config) RESULT(xmap)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: src_idxlist
    TYPE(xt_idxlist), INTENT(in) :: dst_idxlist
    INTEGER, INTENT(in) :: inter_comm, intra_comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_xmap) :: xmap

    INTERFACE
      FUNCTION xt_xmap_dist_dir_intercomm_custom_new_f(src_idxlist, &
           dst_idxlist, inter_comm, intra_comm, config) &
           BIND(C, name='xt_xmap_dist_dir_intercomm_custom_new_f') &
           RESULT(xmap_ptr)
        IMPORT :: xt_idxlist, xt_xmap, xt_mpi_fint_kind, xt_config, c_ptr
        IMPLICIT NONE
        TYPE(xt_idxlist), INTENT(in) :: src_idxlist, dst_idxlist
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: inter_comm, intra_comm
        TYPE(xt_config), INTENT(in) :: config
        TYPE(c_ptr) :: xmap_ptr
      END FUNCTION xt_xmap_dist_dir_intercomm_custom_new_f
    END INTERFACE

    xmap%cptr = xt_xmap_dist_dir_intercomm_custom_new_f(src_idxlist, &
         dst_idxlist, inter_comm, intra_comm, config)
  END FUNCTION xt_xmap_dist_dir_intercomm_custom_new

  FUNCTION xt_xmap_get_num_destinations(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_num_destinations_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_num_destinations')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_num_destinations_c
    END INTERFACE
    num = INT(xt_xmap_get_num_destinations_c(xmap%cptr))
  END FUNCTION xt_xmap_get_num_destinations

  FUNCTION xt_xmap_get_num_sources(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_num_sources_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_num_sources')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_num_sources_c
    END INTERFACE
    num = INT(xt_xmap_get_num_sources_c(xmap%cptr))
  END FUNCTION xt_xmap_get_num_sources

  SUBROUTINE xt_xmap_get_destination_ranks(xmap, ranks)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(out) :: ranks(*)
    INTERFACE
      SUBROUTINE xt_xmap_get_destination_ranks_c(xmap, ranks) &
           BIND(c, name='xt_xmap_get_destination_ranks')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int), INTENT(out) :: ranks(*)
      END SUBROUTINE xt_xmap_get_destination_ranks_c
    END INTERFACE
    CALL xt_xmap_get_destination_ranks_c(xmap%cptr, ranks)
  END SUBROUTINE xt_xmap_get_destination_ranks

  SUBROUTINE xt_xmap_get_source_ranks(xmap, ranks)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(out) :: ranks(*)
    INTERFACE
      SUBROUTINE xt_xmap_get_source_ranks_c(xmap, ranks) &
           BIND(c, name='xt_xmap_get_source_ranks')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int), INTENT(out) :: ranks(*)
      END SUBROUTINE xt_xmap_get_source_ranks_c
    END INTERFACE
    CALL xt_xmap_get_source_ranks_c(xmap%cptr, ranks)
  END SUBROUTINE xt_xmap_get_source_ranks

  FUNCTION xt_xmap_get_max_src_pos(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_max_src_pos_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_max_src_pos')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_max_src_pos_c
    END INTERFACE
    num = INT(xt_xmap_get_max_src_pos_c(xmap%cptr))
  END FUNCTION xt_xmap_get_max_src_pos

  FUNCTION xt_xmap_get_max_dst_pos(xmap) RESULT(num)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER :: num
    INTERFACE
      FUNCTION xt_xmap_get_max_dst_pos_c(xmap) RESULT(num) &
           BIND(c, name='xt_xmap_get_max_dst_pos')
        IMPORT :: c_ptr, c_int
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(c_int) :: num
      END FUNCTION xt_xmap_get_max_dst_pos_c
    END INTERFACE
    num = INT(xt_xmap_get_max_dst_pos_c(xmap%cptr))
  END FUNCTION xt_xmap_get_max_dst_pos

  FUNCTION xt_xmap_reorder(xmap, reorder_type) RESULT(xmap_reorder)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(xt_reorder_type_kind), INTENT(in) :: reorder_type
    TYPE(xt_xmap) :: xmap_reorder
    ! fixme: evil hack because gfortran with all warnings enabled throws a
    ! hissy fit for dummy arguments of type integer(xt_reorder_kind)
    ! of bind(c) routines, but the value of xt_reorder_kind always
    ! matches c_int for that compiler anyway
#ifdef __GNUC__
    LOGICAL, PARAMETER :: static_assert_xt_reorder_type_kind &
         = xt_reorder_type_kind == c_int
    INTEGER, PARAMETER :: assert_check &
         = 1 / MERGE(1, 0, static_assert_xt_reorder_type_kind)
#  define xt_reorder_type_kind c_int
#endif
    INTERFACE
      FUNCTION xt_xmap_reorder_c(xmap, reorder_type) &
           BIND(C, name='xt_xmap_reorder') RESULT(xmap_reorder_ptr)
        IMPORT:: xt_reorder_type_kind, c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        INTEGER(xt_reorder_type_kind), VALUE, INTENT(in) :: reorder_type
        TYPE(c_ptr) :: xmap_reorder_ptr
      END FUNCTION xt_xmap_reorder_c
    END INTERFACE
#ifdef __GNUC__
#  define UNUSED(x) IF (SIZE( (/(x)/) ) < 0) CONTINUE
    UNUSED(assert_check)
#endif

    IF (reorder_type < 0_xt_reorder_type_kind .OR. &
         reorder_type > HUGE(1_c_int)) &
         CALL xt_abort("invalid reorder type", filename, __LINE__)
    xmap_reorder%cptr = xt_xmap_reorder_c(xmap%cptr, reorder_type)
  ! fixme: undo effect of above hack
#ifdef __GNUC__
#  undef xt_reorder_type_kind
#endif
  END FUNCTION xt_xmap_reorder

  FUNCTION xt_xmap_update_positions(xmap, src_positions, dst_positions) &
      RESULT(xmap_updated)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, TARGET, INTENT(in) :: src_positions(*)
    INTEGER, TARGET, INTENT(in) :: dst_positions(*)
    TYPE(xt_xmap) :: xmap_updated
    INTEGER(c_int), TARGET, ALLOCATABLE :: src_positions_c(:), dst_positions_c(:)
    TYPE(c_ptr) :: src_positions_p, dst_positions_p
    INTERFACE
      FUNCTION xt_xmap_update_positions_c(xmap, src_positions, dst_positions) &
           BIND(C, name='xt_xmap_update_positions') RESULT(xmap_updated_ptr)
        IMPORT:: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap, src_positions, dst_positions
        TYPE(c_ptr) :: xmap_updated_ptr
      END FUNCTION xt_xmap_update_positions_c
    END INTERFACE

    IF (c_int == KIND(1)) THEN
      src_positions_p = C_LOC(src_positions)
      dst_positions_p = C_LOC(dst_positions)
    ELSE
      CALL arg2ci(xt_xmap_get_max_src_pos(xmap), src_positions, src_positions_c)
      src_positions_p = C_LOC(src_positions_c)
      CALL arg2ci(xt_xmap_get_max_dst_pos(xmap), dst_positions, dst_positions_c)
      dst_positions_p = C_LOC(dst_positions_c)
    END IF
    xmap_updated%cptr = &
      xt_xmap_update_positions_c(xmap%cptr, src_positions_p, dst_positions_p)
    CONTAINS
      SUBROUTINE arg2ci(n, arg, argc)
        INTEGER, INTENT(in) :: n, arg(*)
        INTEGER(c_int), ALLOCATABLE, INTENT(inout) :: argc(:)
        INTEGER :: i
        ALLOCATE(argc(n))
        DO i = 1, n
          argc(i) = INT(arg(i), c_int)
        END DO
      END SUBROUTINE arg2ci
  END FUNCTION xt_xmap_update_positions

  FUNCTION xt_xmap_spread_a1d(xmap, src_displacements, dst_displacements) &
      RESULT(xmap_spread)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: src_displacements(:)
    INTEGER, INTENT(in) :: dst_displacements(:)
    TYPE(xt_xmap) :: xmap_spread
    INTEGER :: num_repetitions
    INTEGER(i8) :: num_repetitions_i8
    num_repetitions = SIZE(src_displacements)
    IF (num_repetitions /= SIZE(dst_displacements)) &
         CALL xt_abort("invalid number of repetitions", filename, __LINE__)
    num_repetitions_i8 = INT(num_repetitions, i8)
    xmap_spread = &
      xt_xmap_spread( &
        xmap, num_repetitions_i8, src_displacements, dst_displacements);
  END FUNCTION xt_xmap_spread_a1d

  FUNCTION xt_xmap_spread_i2_a1d(xmap, num_repetitions, src_displacements, &
                                 dst_displacements) &
      RESULT(xmap_spread)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i2), INTENT(in) :: num_repetitions
    INTEGER, INTENT(in) :: src_displacements(num_repetitions)
    INTEGER, INTENT(in) :: dst_displacements(num_repetitions)
    TYPE(xt_xmap) :: xmap_spread
    INTEGER(i8) :: num_repetitions_i8
    num_repetitions_i8 = INT(num_repetitions, i8)
    xmap_spread = &
      xt_xmap_spread( &
        xmap, num_repetitions_i8, src_displacements, dst_displacements);
  END FUNCTION xt_xmap_spread_i2_a1d

  FUNCTION xt_xmap_spread_i4_a1d(xmap, num_repetitions, src_displacements, &
                                 dst_displacements) &
      RESULT(xmap_spread)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i4), INTENT(in) :: num_repetitions
    INTEGER, INTENT(in) :: src_displacements(num_repetitions)
    INTEGER, INTENT(in) :: dst_displacements(num_repetitions)
    TYPE(xt_xmap) :: xmap_spread
    INTEGER(i8) :: num_repetitions_i8
    num_repetitions_i8 = INT(num_repetitions, i8)
    xmap_spread = &
      xt_xmap_spread( &
        xmap, num_repetitions_i8, src_displacements, dst_displacements);
  END FUNCTION xt_xmap_spread_i4_a1d

  FUNCTION xt_xmap_spread_i8_a1d(xmap, num_repetitions, src_displacements, &
                                 dst_displacements) &
      RESULT(xmap_spread)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i8), INTENT(in) :: num_repetitions
    INTEGER, TARGET, INTENT(in) :: src_displacements(num_repetitions)
    INTEGER, TARGET, INTENT(in) :: dst_displacements(num_repetitions)
    INTEGER(c_int) :: num_repetitions_c
    TYPE(xt_xmap) :: xmap_spread
    INTEGER(c_int), TARGET, ALLOCATABLE :: &
      src_displacements_c(:), dst_displacements_c(:)
    TYPE(c_ptr) :: src_displacements_p, dst_displacements_p
    INTERFACE
      FUNCTION xt_xmap_spread_c(xmap, num_repetitions, src_displacements, &
                                dst_displacements) &
           BIND(C, name='xt_xmap_spread') RESULT(xmap_spread_ptr)
        IMPORT:: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: &
          xmap, src_displacements, dst_displacements
        INTEGER(c_int), VALUE, INTENT(in) :: num_repetitions
        TYPE(c_ptr) :: xmap_spread_ptr
      END FUNCTION xt_xmap_spread_c
    END INTERFACE
    IF (num_repetitions < 0_c_int .OR. &
        num_repetitions > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_repetitions_c = INT(num_repetitions, c_int)
    IF (c_int == KIND(1)) THEN
      src_displacements_p = C_LOC(src_displacements)
      dst_displacements_p = C_LOC(dst_displacements)
    ELSE
      CALL arg2ci(src_displacements, src_displacements_c)
      src_displacements_p = C_LOC(src_displacements_c)
      CALL arg2ci(dst_displacements, dst_displacements_c)
      dst_displacements_p = C_LOC(dst_displacements_c)
    END IF
    xmap_spread%cptr = &
      xt_xmap_spread_c( &
        xmap%cptr, num_repetitions_c, src_displacements_p, dst_displacements_p)
  CONTAINS
    SUBROUTINE arg2ci(arg, argc)
      INTEGER, INTENT(in) :: arg(*)
      INTEGER(c_int), ALLOCATABLE, INTENT(inout) :: argc(:)
      INTEGER :: i, n
      n = INT(num_repetitions)
      ALLOCATE(argc(n))
      DO i = 1, n
        argc(i) = INT(arg(i), c_int)
      END DO
    END SUBROUTINE arg2ci
  END FUNCTION xt_xmap_spread_i8_a1d

  FUNCTION xt_xmap_get_out_iterator(xmap) RESULT(iter)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_xmap_iter) :: iter
    INTERFACE
      FUNCTION xt_xmap_get_out_iterator_c(xmap) RESULT(cptr) &
           BIND(c, name='xt_xmap_get_out_iterator')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        TYPE(c_ptr) :: cptr
      END FUNCTION xt_xmap_get_out_iterator_c
    END INTERFACE
    iter%cptr = xt_xmap_get_out_iterator_c(xmap%cptr)
  END FUNCTION xt_xmap_get_out_iterator

  FUNCTION xt_xmap_get_in_iterator(xmap) RESULT(iter)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_xmap_iter) :: iter
    INTERFACE
      FUNCTION xt_xmap_get_in_iterator_c(xmap) RESULT(cptr) &
           BIND(c, name='xt_xmap_get_in_iterator')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: xmap
        TYPE(c_ptr) :: cptr
      END FUNCTION xt_xmap_get_in_iterator_c
    END INTERFACE
    iter%cptr = xt_xmap_get_in_iterator_c(xmap%cptr)
  END FUNCTION xt_xmap_get_in_iterator

  FUNCTION xt_xmap_iterator_is_null(iter) RESULT(p)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(iter%cptr)
  END FUNCTION xt_xmap_iterator_is_null

  FUNCTION xt_xmap_iterator_next(iter) RESULT(avail)
    TYPE(xt_xmap_iter), INTENT(inout) :: iter
    LOGICAL :: avail
    INTERFACE
      FUNCTION xt_xmap_iterator_next_c(iter) RESULT(avail) &
           BIND(c, name='xt_xmap_iterator_next')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: avail
      END FUNCTION xt_xmap_iterator_next_c
    END INTERFACE
    avail = xt_xmap_iterator_next_c(iter%cptr) /= 0
  END FUNCTION xt_xmap_iterator_next

  FUNCTION xt_xmap_iterator_get_rank(iter) RESULT(rank)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER :: rank
    INTERFACE
      FUNCTION xt_xmap_iterator_get_rank_c(iter) RESULT(rank) &
           BIND(c, name='xt_xmap_iterator_get_rank')
        IMPORT :: c_ptr, c_int
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        INTEGER(c_int) :: rank
      END FUNCTION xt_xmap_iterator_get_rank_c
    END INTERFACE
    rank = INT(xt_xmap_iterator_get_rank_c(iter%cptr))
  END FUNCTION xt_xmap_iterator_get_rank

  !> note: result is read-only
  FUNCTION xt_xmap_iterator_get_transfer_pos(iter) RESULT(transfer_pos)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER(c_int), POINTER :: transfer_pos(:)

    INTERFACE
      FUNCTION xt_xmap_iterator_get_transfer_pos_c(iter) RESULT(transfer_pos) &
           BIND(c, name='xt_xmap_iterator_get_transfer_pos')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        TYPE(c_ptr) :: transfer_pos
      END FUNCTION xt_xmap_iterator_get_transfer_pos_c
    END INTERFACE
    INTEGER :: n(1)
    TYPE(c_ptr) :: transfer_pos_cptr
    NULLIFY(transfer_pos)
    n(1) = INT(xt_xmap_iterator_get_num_transfer_pos_c(iter%cptr))
    transfer_pos_cptr = xt_xmap_iterator_get_transfer_pos_c(iter%cptr)
    CALL C_F_POINTER(transfer_pos_cptr, transfer_pos, n)
  END FUNCTION xt_xmap_iterator_get_transfer_pos

  FUNCTION xt_xmap_iterator_get_num_transfer_pos(iter) RESULT(num)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER :: num
    num = INT(xt_xmap_iterator_get_num_transfer_pos_c(iter%cptr))
  END FUNCTION xt_xmap_iterator_get_num_transfer_pos

  !> note: result is read-only
  FUNCTION xt_xmap_iterator_get_transfer_pos_ext(iter) RESULT(transfer_pos_ext)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    TYPE(xt_pos_ext), POINTER :: transfer_pos_ext(:)

    INTERFACE
      FUNCTION xt_xmap_iterator_get_transfer_pos_ext_c(iter) &
           RESULT(transfer_pos_ext) &
           BIND(c, name='xt_xmap_iterator_get_transfer_pos_ext')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
        TYPE(c_ptr) :: transfer_pos_ext
      END FUNCTION xt_xmap_iterator_get_transfer_pos_ext_c
    END INTERFACE
    INTEGER :: n(1)
    TYPE(c_ptr) :: transfer_pos_ext_cptr
    NULLIFY(transfer_pos_ext)
    n(1) = INT(xt_xmap_iterator_get_num_transfer_pos_ext_c(iter%cptr))
    transfer_pos_ext_cptr = xt_xmap_iterator_get_transfer_pos_ext_c(iter%cptr)
    CALL C_F_POINTER(transfer_pos_ext_cptr, transfer_pos_ext, n)
  END FUNCTION xt_xmap_iterator_get_transfer_pos_ext

  FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext(iter) RESULT(num)
    TYPE(xt_xmap_iter), INTENT(in) :: iter
    INTEGER :: num
    num = INT(xt_xmap_iterator_get_num_transfer_pos_ext_c(iter%cptr))
  END FUNCTION xt_xmap_iterator_get_num_transfer_pos_ext

  SUBROUTINE xt_xmap_iterator_delete(iter)
    TYPE(xt_xmap_iter), INTENT(inout) :: iter
    INTERFACE
      SUBROUTINE xt_xmap_iterator_delete_c(iter) &
           BIND(c, name='xt_xmap_iterator_delete')
        IMPORT :: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: iter
      END SUBROUTINE xt_xmap_iterator_delete_c
    END INTERFACE
    CALL xt_xmap_iterator_delete_c(iter%cptr)
    iter%cptr = c_null_ptr
  END SUBROUTINE xt_xmap_iterator_delete
END MODULE xt_xmap_abstract

MODULE xt_xmap_rename
  USE xt_xmap_abstract, ONLY: xt_xmap_all2all_orig_new => xt_xmap_all2all_new, &
       xt_xmap_all2all_custom_new, &
       xt_xmap_dist_dir_orig_new => xt_xmap_dist_dir_new, &
       xt_xmap_dist_dir_custom_new, &
       xt_xmap_dist_dir_intercomm_orig_new => xt_xmap_dist_dir_intercomm_new, &
       xt_xmap_dist_dir_intercomm_custom_new
  IMPLICIT NONE
  PRIVATE
  INTERFACE xt_xmap_all2all_new
    MODULE PROCEDURE xt_xmap_all2all_orig_new
    MODULE PROCEDURE xt_xmap_all2all_custom_new
  END INTERFACE xt_xmap_all2all_new
  PUBLIC :: xt_xmap_all2all_new
  INTERFACE xt_xmap_dist_dir_new
    MODULE PROCEDURE xt_xmap_dist_dir_orig_new
    MODULE PROCEDURE xt_xmap_dist_dir_custom_new
  END INTERFACE xt_xmap_dist_dir_new
  PUBLIC :: xt_xmap_dist_dir_new
  INTERFACE xt_xmap_dist_dir_intercomm_new
    MODULE PROCEDURE xt_xmap_dist_dir_intercomm_orig_new
    MODULE PROCEDURE xt_xmap_dist_dir_intercomm_custom_new
  END INTERFACE xt_xmap_dist_dir_intercomm_new
  PUBLIC :: xt_xmap_dist_dir_intercomm_new
END MODULE xt_xmap_rename
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
