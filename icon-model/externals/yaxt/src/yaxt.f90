!>
!! @file yaxt.f90
!! @brief Fortran interface to yaxt implementation
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

MODULE yaxt
  !
  ! Fortran interface to yaxt implementation
  !
  USE xt_core, ONLY: i4, xt_int_kind, xt_int_mpidt, &
       xt_abort, xt_initialize, xt_get_default_comm, xt_stripe, xt_bounds, &
       char, xt_finalize, xt_initialized, xt_finalized, xt_slice_c_loc, &
       xt_pos_ext, OPERATOR(/=), OPERATOR(==), &
       xt_set_abort_handler => set_abort_handler, xt_restore_default_abort_hndl
  USE xt_config_f, ONLY: xt_config, xt_config_new, xt_config_delete, &
       xt_exchanger_id_by_name, &
       xt_config_get_exchange_method, xt_config_set_exchange_method, &
       xt_exchanger_irecv_send, xt_exchanger_irecv_isend, &
       xt_exchanger_irecv_isend_packed, xt_exchanger_mix_isend_irecv, &
       xt_exchanger_neigh_alltoall, &
       xt_config_get_idxvec_autoconvert_size, &
       xt_config_set_idxvec_autoconvert_size, &
       xt_config_get_redist_mthread_mode, &
       xt_config_set_redist_mthread_mode, &
       xt_exchanger_irecv_isend_ddt_packed
  USE xt_sort, ONLY: xt_sort_int, xt_sort_index, xt_sort_idxpos, &
       xt_sort_permutation, xt_assign_id_map
  USE xt_idxlist_abstract, ONLY: &
       xt_idxlist_f2c, xt_idxlist_c2f, xt_idxlist, &
       xt_idxlist_delete, xt_idxlist_get_pack_size, &
       xt_idxlist_pack => xt_idxlist_pack_f, &
       xt_idxlist_unpack => xt_idxlist_unpack_f, xt_idxlist_copy, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_idxlist_get_positions_of_indices, xt_idxlist_get_index_stripes, &
       xt_idxlist_get_bounding_box, xt_idxlist_get_intersection, &
       xt_idxlist_get_num_indices, xt_is_null, &
       xt_idxlist_get_indices, xt_idxlist_get_indices_const, &
       xt_idxlist_get_pos_exts_of_index_stripes
  USE xt_idxvec, ONLY: xt_idxvec_from_stripes_new, xt_idxvec_new
  USE xt_idxstripes, ONLY: xt_idxstripes_new, xt_idxstripes_from_idxlist_new
  USE xt_idxsection, ONLY: xt_idxsection_new, xt_idxfsection_new
  USE xt_idxlist_collection, ONLY: xt_idxlist_collection_new
  USE xt_xmap_abstract, ONLY: xt_xmap, xt_xmap_c2f, xt_xmap_f2c, &
       xt_xmap_delete, xt_xmap_get_num_destinations, xt_xmap_get_num_sources, &
       xt_xmap_copy, xt_xmap_get_destination_ranks, xt_xmap_get_source_ranks, &
       xt_xmap_get_max_src_pos, xt_xmap_get_max_dst_pos, &
       xt_xmap_all2all_custom_new, xt_xmap_dist_dir_custom_new, xt_is_null, &
       xt_xmap_dist_dir_intercomm_custom_new, &
       xt_xmap_iter, xt_xmap_get_in_iterator, xt_xmap_get_out_iterator, &
       xt_xmap_iterator_next, xt_xmap_iterator_get_rank, &
       xt_xmap_iterator_get_transfer_pos, &
       xt_xmap_iterator_get_num_transfer_pos, xt_xmap_iterator_delete, &
       xt_xmap_iterator_get_transfer_pos_ext, &
       xt_xmap_iterator_get_num_transfer_pos_ext, xt_xmap_reorder, &
       xt_reorder_type_kind, XT_REORDER_NONE, XT_REORDER_SEND_UP, &
       XT_REORDER_RECV_UP, xt_xmap_update_positions, xt_xmap_spread
  USE xt_xmap_rename, ONLY: xt_xmap_all2all_new, xt_xmap_dist_dir_new, &
       xt_xmap_dist_dir_intercomm_new
  USE xt_xmap_intersection, ONLY: xt_xmap_intersection_new, &
       xt_xmap_intersection_ext_new, xt_com_list, &
       xt_xmap_intersection_pos_new, xt_com_pos
  USE xt_redist_base, ONLY: xt_redist, xt_redist_c2f, xt_redist_f2c, &
       xt_redist_copy, xt_redist_p2p_custom_new, &
       xt_redist_delete, xt_redist_s_exchange1, xt_redist_s_exchange, &
       xt_redist_p2p_off_new, xt_redist_p2p_off_custom_new, &
       xt_redist_p2p_blocks_new, xt_redist_p2p_blocks_custom_new, &
       xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_off_custom_new, &
       xt_redist_collection_static_new, xt_redist_collection_new, &
       xt_redist_repeat_new, xt_is_null, xt_redist_get_mpi_comm, &
       xt_offset_ext, xt_redist_p2p_ext_new, xt_redist_msg, &
       xt_aoffset_ext, xt_redist_p2p_aext_new, &
       xt_redist_a_exchange1, xt_redist_a_exchange, &
       xt_redist_single_array_base_new, &
       xt_redist_single_array_base_custom_new, &
       xt_redist_get_recv_mpi_datatype, xt_redist_get_send_mpi_datatype, &
       xt_redist_get_num_send_msg, xt_redist_get_num_recv_msg
  USE xt_redist_rename, ONLY: xt_redist_p2p_new
  USE xt_redist_int_i2, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_int_i4, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_int_i8, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_logical, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_real_sp, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_redist_real_dp, ONLY: xt_redist_s_exchange, xt_redist_a_exchange
  USE xt_requests, ONLY: xt_request, xt_request_wait, xt_request_test, &
       xt_is_null, xt_request_null
  USE xt_mpi, ONLY: xt_mpi_comm_mark_exclusive

  USE iso_c_binding, ONLY: c_int, c_ptr, c_null_ptr, c_loc
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: xt_initialize, xt_finalize, xt_abort, xt_get_default_comm, &
       xt_initialized, xt_finalized, &
       xt_set_abort_handler, xt_restore_default_abort_hndl, &
       xt_int_kind, xt_int_mpidt, xt_stripe, xt_bounds, xt_pos_ext, &
       xt_idxlist, xt_idxlist_delete,  &
       xt_idxlist_c2f, xt_idxlist_f2c, xt_idxlist_copy, &
       xt_idxlist_get_pack_size, xt_idxlist_pack, xt_idxlist_unpack, &
       xt_idxlist_get_intersection, xt_idxlist_get_positions_of_indices, &
       xt_idxlist_get_index_at_position, xt_idxlist_get_indices_at_positions, &
       xt_idxlist_get_position_of_index, xt_idxlist_get_position_of_index_off, &
       xt_idxlist_get_num_indices, &
       xt_idxlist_get_indices, xt_idxlist_get_indices_const, &
       xt_idxlist_get_index_stripes, xt_idxlist_get_bounding_box, &
       xt_idxlist_get_pos_exts_of_index_stripes, &
       xt_idxvec_new, xt_idxvec_from_stripes_new, &
       xt_idxstripes_new, xt_idxstripes_from_idxlist_new, &
       xt_idxlist_collection_new, xt_idxsection_new, xt_idxfsection_new, &
       xt_idxmod_new, xt_xmap, xt_xmap_c2f, xt_xmap_f2c, xt_xmap_copy, &
       xt_xmap_delete, xt_xmap_get_num_destinations, &
       xt_xmap_get_num_sources, xt_xmap_get_destination_ranks, &
       xt_xmap_get_source_ranks, &
       xt_xmap_get_max_src_pos, xt_xmap_get_max_dst_pos, &
       xt_xmap_intersection_new, xt_xmap_intersection_ext_new, &
       xt_xmap_intersection_pos_new, &
       xt_xmap_reorder, xt_reorder_type_kind, &
       XT_REORDER_NONE, XT_REORDER_SEND_UP, XT_REORDER_RECV_UP, &
       xt_xmap_update_positions, xt_xmap_spread, &
       xt_com_list, xt_com_pos, &
       xt_xmap_iter, xt_xmap_get_in_iterator, xt_xmap_get_out_iterator, &
       xt_xmap_iterator_next, xt_xmap_iterator_get_rank, &
       xt_xmap_iterator_get_transfer_pos, &
       xt_xmap_iterator_get_num_transfer_pos, xt_xmap_iterator_delete, &
       xt_xmap_iterator_get_transfer_pos_ext, &
       xt_xmap_iterator_get_num_transfer_pos_ext, &
       xt_redist, xt_redist_f2c, xt_redist_c2f, &
       xt_redist_msg, &
       xt_offset_ext, xt_redist_p2p_ext_new, &
       xt_aoffset_ext, xt_redist_p2p_aext_new, &
       xt_redist_repeat_new, xt_redist_copy, xt_redist_delete, &
       xt_redist_s_exchange1, xt_redist_s_exchange, &
       xt_redist_get_mpi_comm, &
       xt_redist_get_recv_mpi_datatype, xt_redist_get_send_mpi_datatype, &
       xt_redist_get_num_send_msg, xt_redist_get_num_recv_msg, &
       xt_idxempty_new, xt_redist_collection_static_new, &
       xt_redist_collection_new, xt_slice_c_loc, &
       xt_mpi_comm_mark_exclusive, &
       xt_is_null, &
       OPERATOR(==), OPERATOR(/=), char, &
       xt_request, xt_request_null, &
       xt_request_wait, xt_request_test, xt_redist_a_exchange1, &
       xt_redist_a_exchange

  PUBLIC :: xt_sort_int, xt_sort_index, xt_sort_idxpos, xt_sort_permutation, &
       xt_assign_id_map
  PUBLIC :: xt_xmap_all2all_new, xt_xmap_all2all_custom_new
  PUBLIC :: xt_xmap_dist_dir_new, xt_xmap_dist_dir_custom_new
  PUBLIC :: xt_xmap_dist_dir_intercomm_new, &
       xt_xmap_dist_dir_intercomm_custom_new

  PUBLIC :: xt_config, xt_config_new, xt_config_delete, &
       xt_exchanger_id_by_name, &
       xt_config_get_exchange_method, xt_config_set_exchange_method, &
       xt_exchanger_irecv_send, xt_exchanger_irecv_isend, &
       xt_exchanger_irecv_isend_packed, xt_exchanger_mix_isend_irecv, &
       xt_exchanger_neigh_alltoall, &
       xt_config_get_idxvec_autoconvert_size, &
       xt_config_set_idxvec_autoconvert_size, &
       xt_config_get_redist_mthread_mode, &
       xt_config_set_redist_mthread_mode, &
       xt_exchanger_irecv_isend_ddt_packed

  PUBLIC :: xt_redist_p2p_new, xt_redist_p2p_custom_new
  PUBLIC :: xt_redist_p2p_off_new, xt_redist_p2p_off_custom_new
  PUBLIC :: xt_redist_p2p_blocks_new, xt_redist_p2p_blocks_custom_new
  PUBLIC :: xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_off_custom_new
  PUBLIC :: xt_redist_single_array_base_new, &
       xt_redist_single_array_base_custom_new

  INTERFACE OPERATOR(==)
    MODULE PROCEDURE xt_bounds_eq
  END INTERFACE OPERATOR(==)

  INTERFACE OPERATOR(/=)
    MODULE PROCEDURE xt_bounds_ne
  END INTERFACE OPERATOR(/=)

  TYPE, BIND(C), PUBLIC :: xt_modifier
    TYPE(xt_idxlist) :: extract
    TYPE(xt_idxlist) :: subst
    INTEGER(c_int) :: mask
  END TYPE xt_modifier

  INTERFACE
    FUNCTION xt_idxmod_new_c(patch, modifier, modifier_num, mstate_ptr) &
         BIND(C, name='xt_idxmod_new') RESULT(res)
      IMPORT :: xt_int_kind, xt_modifier, c_ptr, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: patch
      INTEGER(c_int), VALUE, INTENT(in) :: modifier_num
      TYPE(xt_modifier), INTENT(in) :: modifier(modifier_num)
      TYPE(c_ptr), VALUE, INTENT(in) :: mstate_ptr
      TYPE(c_ptr) :: res
    END FUNCTION xt_idxmod_new_c
  END INTERFACE

  INTERFACE xt_idxmod_new
    MODULE PROCEDURE xt_idxmod_new_a1d
    MODULE PROCEDURE xt_idxmod_new_a1d_a1d
    MODULE PROCEDURE xt_idxmod_new_a1d_i4
    MODULE PROCEDURE xt_idxmod_new_a1d_i4_a1d
    MODULE PROCEDURE xt_idxmod_new_a1d_i4_a2d
  END INTERFACE xt_idxmod_new
  CHARACTER(len=*), PARAMETER :: filename = 'yaxt.f90'
CONTAINS

  FUNCTION xt_idxempty_new() RESULT(res)
    IMPLICIT NONE
    TYPE(Xt_idxlist) :: res

    INTERFACE
      FUNCTION xt_idxempty_new_c() &
           BIND(C, name='xt_idxempty_new') RESULT(res_ptr)
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_idxempty_new_c
    END INTERFACE

    res = xt_idxlist_c2f(xt_idxempty_new_c())

  END FUNCTION xt_idxempty_new

  ELEMENTAL FUNCTION xt_bounds_eq(a, b) RESULT(a_equals_b)
    TYPE(xt_bounds), INTENT(in) :: a, b
    LOGICAL :: a_equals_b
    a_equals_b = a%size == b%size .AND. a%start == b%start
  END FUNCTION xt_bounds_eq

  ELEMENTAL FUNCTION xt_bounds_ne(a, b) RESULT(a_equals_b)
    TYPE(xt_bounds), INTENT(in) :: a, b
    LOGICAL :: a_equals_b
    a_equals_b = a%size /= b%size .OR. a%start /= b%start
  END FUNCTION xt_bounds_ne

  FUNCTION xt_idxmod_new_a1d(patch, modifier) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(:)
    TYPE(Xt_idxlist) :: res

    INTEGER :: num_modifier
    INTEGER(c_int) :: num_modifier_c
    num_modifier = SIZE(modifier)
    IF (num_modifier > HUGE(1_c_int)) &
         CALL xt_abort("number of modifiers too high", filename, __LINE__)
    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, c_null_ptr))
  END FUNCTION xt_idxmod_new_a1d

  FUNCTION xt_idxmod_new_a1d_a1d(patch, modifier, mstate) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(:)
    INTEGER(c_int), TARGET, INTENT(inout) :: mstate(*)
    TYPE(Xt_idxlist) :: res

    INTEGER :: num_modifier
    INTEGER(c_int) :: num_modifier_c
    num_modifier = SIZE(modifier)
    IF (num_modifier > HUGE(1_c_int)) &
         CALL xt_abort("number of modifiers too high", filename, __LINE__)
    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, C_LOC(mstate)))
  END FUNCTION xt_idxmod_new_a1d_a1d

  FUNCTION xt_idxmod_new_a1d_i4(patch, modifier, num_modifier) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(*)
    INTEGER(i4), INTENT(in) :: num_modifier
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_modifier_c

    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, c_null_ptr))
  END FUNCTION xt_idxmod_new_a1d_i4

  FUNCTION xt_idxmod_new_a1d_i4_a1d(patch, modifier, num_modifier, mstate) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(*)
    INTEGER(i4), INTENT(in) :: num_modifier
    INTEGER(c_int), TARGET, INTENT(inout) :: mstate(*)
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_modifier_c

    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, C_LOC(mstate)))
  END FUNCTION xt_idxmod_new_a1d_i4_a1d

  FUNCTION xt_idxmod_new_a1d_i4_a2d(patch, modifier, num_modifier, mstate) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_idxlist), INTENT(in) :: patch
    TYPE(xt_modifier), INTENT(in) :: modifier(*)
    INTEGER(i4), INTENT(in) :: num_modifier
    INTEGER(c_int), TARGET, INTENT(inout) :: mstate(1,*)
    TYPE(Xt_idxlist) :: res
    INTEGER(c_int) :: num_modifier_c

    num_modifier_c = INT(num_modifier, c_int)
    res = xt_idxlist_c2f(xt_idxmod_new_c(xt_idxlist_f2c(patch), modifier, &
         num_modifier_c, C_LOC(mstate)))
  END FUNCTION xt_idxmod_new_a1d_i4_a2d

END MODULE yaxt
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
