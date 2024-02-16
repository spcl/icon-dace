!>
!! @file xt_redist_f.f90
!! @brief xt_redist-related procedures of Fortran interface
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
!! @example test_redist_collection_displace_f.f90
!! @example test_redist_common_f.f90
!! @example test_redist_p2p_f.f90
!! @example test_redist_p2p_parallel_f.f90

#include "fc_feature_defs.inc"
MODULE xt_redist_base
  USE xt_core, ONLY: xt_abort, xt_mpi_fint_kind, i2, i4, i8
  USE xt_config_f, ONLY: xt_config, xt_config_f2c
  USE xt_xmap_abstract, ONLY: xt_xmap
  USE iso_c_binding, ONLY: c_int, c_null_ptr, c_ptr, c_loc, c_associated
  USE xt_mpi, ONLY: mpi_address_kind
  USE xt_requests, ONLY: xt_request
#include "xt_slice_c_loc.inc"
  IMPLICIT NONE
  PRIVATE
  ! note: this type must not be extended to contain any other
  ! components, its memory pattern has to match void * exactly, which
  ! it does because of C constraints
  TYPE, BIND(C), PUBLIC :: xt_redist
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_redist

  TYPE, BIND(c), PUBLIC :: xt_offset_ext
    INTEGER(c_int) :: start, size, stride
  END TYPE xt_offset_ext

  TYPE, BIND(c), PUBLIC :: xt_redist_msg
    INTEGER(xt_mpi_fint_kind) :: rank, datatype
  END TYPE xt_redist_msg

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_redist_f2c(redist) BIND(c, name='xt_redist_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_redist
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(in) :: redist
      TYPE(c_ptr) :: p
    END FUNCTION xt_redist_f2c
  END INTERFACE

  INTERFACE xt_redist_delete
    MODULE PROCEDURE xt_redist_delete_1
    MODULE PROCEDURE xt_redist_delete_a1d
  END INTERFACE xt_redist_delete

  INTERFACE
    SUBROUTINE xt_redist_delete_c(redist) &
         BIND(C, name='xt_redist_delete')
      IMPORT :: c_ptr
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: redist
    END SUBROUTINE xt_redist_delete_c
  END INTERFACE

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_redist_is_null
  END INTERFACE xt_is_null

  INTERFACE xt_redist_s_exchange
    MODULE PROCEDURE xt_redist_s_exchange1
    MODULE PROCEDURE xt_redist_s_exchange_a1d
    MODULE PROCEDURE xt_redist_s_exchange_i2_a1d
    MODULE PROCEDURE xt_redist_s_exchange_i4_a1d
    MODULE PROCEDURE xt_redist_s_exchange_i8_a1d
  END INTERFACE xt_redist_s_exchange

  INTERFACE xt_redist_a_exchange
    MODULE PROCEDURE xt_redist_a_exchange1
    MODULE PROCEDURE xt_redist_a_exchange_a1d
    MODULE PROCEDURE xt_redist_a_exchange_i2_a1d
    MODULE PROCEDURE xt_redist_a_exchange_i4_a1d
    MODULE PROCEDURE xt_redist_a_exchange_i8_a1d
  END INTERFACE xt_redist_a_exchange


  INTERFACE
    SUBROUTINE xt_redist_s_exchange_c(redist, num_ptr, src_data_cptr, &
         dst_data_cptr) BIND(C, name='xt_redist_s_exchange')
      IMPORT:: c_ptr, c_int
      TYPE(c_ptr), VALUE, INTENT(in) :: redist
      INTEGER(c_int), VALUE, INTENT(in) :: num_ptr
      TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    END SUBROUTINE xt_redist_s_exchange_c

    SUBROUTINE xt_redist_a_exchange_c(redist, num_ptr, src_data_cptr, &
         dst_data_cptr, request) BIND(C, name='xt_redist_a_exchange')
      IMPORT:: c_ptr, c_int, xt_request
      TYPE(c_ptr), VALUE, INTENT(in) :: redist
      INTEGER(c_int), VALUE, INTENT(in) :: num_ptr
      TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
      TYPE(xt_request), INTENT(out) :: request
    END SUBROUTINE  xt_redist_a_exchange_c

    FUNCTION xt_redist_get_mpi_comm(redist) &
         BIND(c, name='xt_redist_get_mpi_comm_c2f') RESULT(comm)
      IMPORT :: xt_redist, xt_mpi_fint_kind
      TYPE(xt_redist), INTENT(in) :: redist
      INTEGER(xt_mpi_fint_kind) :: comm
    END FUNCTION xt_redist_get_mpi_comm

    FUNCTION xt_redist_get_num_send_msg_c(redist) RESULT(num_send_msg) &
         BIND(c, name='xt_redist_get_num_send_msg')
      IMPORT :: c_ptr, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: redist
      INTEGER(c_int) :: num_send_msg
    END FUNCTION xt_redist_get_num_send_msg_c

    FUNCTION xt_redist_get_num_recv_msg_c(redist) RESULT(num_recv_msg) &
         BIND(c, name='xt_redist_get_num_recv_msg')
      IMPORT :: c_ptr, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in) :: redist
      INTEGER(c_int) :: num_recv_msg
    END FUNCTION xt_redist_get_num_recv_msg_c

    FUNCTION xt_redist_get_recv_mpi_datatype(redist, rank) &
         BIND(c, name='xt_redist_get_recv_MPI_Datatype_c2f') RESULT(dt)
      IMPORT :: xt_redist, xt_mpi_fint_kind
      TYPE(xt_redist), INTENT(in) :: redist
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: rank
      INTEGER(xt_mpi_fint_kind) :: dt
    END FUNCTION xt_redist_get_recv_mpi_datatype

    FUNCTION xt_redist_get_send_mpi_datatype(redist, rank) &
         BIND(c, name='xt_redist_get_send_MPI_Datatype_c2f') RESULT(dt)
      IMPORT :: xt_redist, xt_mpi_fint_kind
      TYPE(xt_redist), INTENT(in) :: redist
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: rank
      INTEGER(xt_mpi_fint_kind) :: dt
    END FUNCTION xt_redist_get_send_mpi_datatype

    FUNCTION xt_redist_collection_static_new_f(redists_f, num_redists, &
         src_displacements, dst_displacements, comm_f) &
         BIND(C, name='xt_redist_collection_static_new_f') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(in) :: redists_f(*)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: num_redists
      INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(*), &
           dst_displacements(*)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm_f
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_collection_static_new_f

    FUNCTION xt_redist_collection_static_custom_new_f(redists_f, num_redists, &
         src_displacements, dst_displacements, comm_f, config) &
         BIND(C, name='xt_redist_collection_static_custom_new_f') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind, xt_config
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(in) :: redists_f(*)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: num_redists
      INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(*), &
           dst_displacements(*)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: comm_f
      TYPE(xt_config), INTENT(in) :: config
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_collection_static_custom_new_f

    FUNCTION xt_redist_collection_new_f(redists_f, num_redists, cache_size, &
         comm_f) BIND(C, name='xt_redist_collection_new_f') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(in) :: redists_f(*)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: &
           num_redists, cache_size, comm_f
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_collection_new_f

    FUNCTION xt_redist_collection_custom_new_f(redists_f, num_redists, &
         cache_size, comm_f, config) &
         BIND(C, name='xt_redist_collection_custom_new_f') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind, xt_config
      IMPLICIT NONE
      TYPE(xt_redist), INTENT(in) :: redists_f(*)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: &
           num_redists, cache_size, comm_f
      TYPE(xt_config), INTENT(in) :: config
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_collection_custom_new_f

    FUNCTION xt_redist_p2p_ext_new_c2f(xmap, num_src_ext, src_extents, &
         num_dst_ext, dst_extents, datatype) &
         BIND(c, name='xt_redist_p2p_ext_new_c2f') RESULT(redist)
      IMPORT :: c_int, c_ptr, xt_offset_ext, xt_mpi_fint_kind, xt_xmap
      TYPE(xt_xmap), INTENT(in) :: xmap
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_ext, num_dst_ext
      TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
           dst_extents(num_dst_ext)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
      TYPE(c_ptr) :: redist
    END FUNCTION xt_redist_p2p_ext_new_c2f

    FUNCTION xt_redist_p2p_ext_custom_new_c2f(xmap, num_src_ext, src_extents, &
         num_dst_ext, dst_extents, datatype, config) &
         BIND(c, name='xt_redist_p2p_ext_custom_new_c2f') RESULT(redist)
      IMPORT :: c_int, c_ptr, xt_offset_ext, xt_mpi_fint_kind, xt_xmap, &
           xt_config
      TYPE(xt_xmap), INTENT(in) :: xmap
      INTEGER(c_int), VALUE, INTENT(in) :: num_src_ext, num_dst_ext
      TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
           dst_extents(num_dst_ext)
      INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
      TYPE(xt_config), INTENT(in) :: config
      TYPE(c_ptr) :: redist
    END FUNCTION xt_redist_p2p_ext_custom_new_c2f

    FUNCTION xt_redist_repeat_new_c(redist_f, src_extent, dst_extent, &
         num_repetitions, displacements) &
         BIND(C, name='xt_redist_repeat_new') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in)  :: redist_f
      INTEGER(mpi_address_kind), VALUE, INTENT(in) :: src_extent, dst_extent
      INTEGER(c_int), VALUE, INTENT(in) :: num_repetitions
      INTEGER(c_int), INTENT(in) :: displacements(*)
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_repeat_new_c

    FUNCTION xt_redist_repeat_custom_new_c(redist_f, src_extent, dst_extent, &
         num_repetitions, displacements, config) &
         BIND(C, name='xt_redist_repeat_custom_new') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in)  :: redist_f
      INTEGER(mpi_address_kind), VALUE, INTENT(in) :: src_extent, dst_extent
      INTEGER(c_int), VALUE, INTENT(in) :: num_repetitions
      INTEGER(c_int), INTENT(in) :: displacements(*)
      TYPE(c_ptr), VALUE, INTENT(in) :: config
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_repeat_custom_new_c

    FUNCTION xt_redist_repeat_asym_new_c(redist_f, src_extent, dst_extent, &
         num_repetitions, src_displacements, dst_displacements) &
         BIND(C, name='xt_redist_repeat_asym_new') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in)  :: redist_f
      INTEGER(mpi_address_kind), VALUE, INTENT(in) :: src_extent, dst_extent
      INTEGER(c_int), VALUE, INTENT(in) :: num_repetitions
      INTEGER(c_int), INTENT(in) :: src_displacements(*), dst_displacements(*)
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_repeat_asym_new_c

    FUNCTION xt_redist_repeat_asym_custom_new_c(redist_f, src_extent, &
         dst_extent, num_repetitions, src_displacements, dst_displacements, &
         config) &
         BIND(C, name='xt_redist_repeat_asym_custom_new') RESULT(res)
      IMPORT :: xt_redist, mpi_address_kind, c_ptr, xt_mpi_fint_kind, c_int
      IMPLICIT NONE
      TYPE(c_ptr), VALUE, INTENT(in)  :: redist_f
      INTEGER(mpi_address_kind), VALUE, INTENT(in) :: src_extent, dst_extent
      INTEGER(c_int), VALUE, INTENT(in) :: num_repetitions
      INTEGER(c_int), INTENT(in) :: src_displacements(*), dst_displacements(*)
      TYPE(c_ptr), VALUE, INTENT(in) :: config
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_repeat_asym_custom_new_c

    FUNCTION xt_redist_single_array_base_new_c2f(nsend, nrecv, send_msgs_f, &
         recv_msgs_f, comm_f) &
         BIND(C, name='xt_redist_single_array_base_new_c2f') RESULT(res)
      IMPORT :: c_ptr, c_int, xt_mpi_fint_kind
      IMPLICIT NONE
      INTEGER(c_int), VALUE :: nsend, nrecv
      TYPE(c_ptr), VALUE, INTENT(in) :: send_msgs_f, recv_msgs_f
      INTEGER(xt_mpi_fint_kind), VALUE :: comm_f
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_single_array_base_new_c2f

    FUNCTION xt_redist_single_array_base_custom_new_c2f(nsend, nrecv, &
         send_msgs_f, recv_msgs_f, comm_f, config) &
         BIND(C, name='xt_redist_single_array_base_custom_new_c2f') RESULT(res)
      IMPORT :: c_ptr, c_int, xt_mpi_fint_kind, xt_config
      IMPLICIT NONE
      INTEGER(c_int), VALUE :: nsend, nrecv
      TYPE(c_ptr), VALUE, INTENT(in) :: send_msgs_f, recv_msgs_f
      INTEGER(xt_mpi_fint_kind), VALUE :: comm_f
      TYPE(xt_config), INTENT(in) :: config
      TYPE(c_ptr) :: res
    END FUNCTION xt_redist_single_array_base_custom_new_c2f

  END INTERFACE

  INTERFACE xt_redist_collection_static_new
    MODULE PROCEDURE xt_redist_collection_static_new_a_i_2ak_i
    MODULE PROCEDURE xt_redist_collection_static_new_a_i_2ak_i_cfg
    MODULE PROCEDURE xt_redist_collection_static_new_a_2ak_i
    MODULE PROCEDURE xt_redist_collection_static_new_a_2ak_i_cfg
  END INTERFACE xt_redist_collection_static_new

  INTERFACE xt_redist_collection_static_custom_new
    MODULE PROCEDURE xt_redist_collection_static_new_a_i_2ak_i_cfg
    MODULE PROCEDURE xt_redist_collection_static_new_a_2ak_i_cfg
  END INTERFACE xt_redist_collection_static_custom_new

  INTERFACE xt_redist_collection_new
    MODULE PROCEDURE xt_redist_collection_new_a_3i
    MODULE PROCEDURE xt_redist_collection_new_a_3i_cfg
    MODULE PROCEDURE xt_redist_collection_new_a_2i
    MODULE PROCEDURE xt_redist_collection_new_a_2i_cfg
    MODULE PROCEDURE xt_redist_collection_new_a_i
    MODULE PROCEDURE xt_redist_collection_new_a_i_cfg
  END INTERFACE xt_redist_collection_new

  INTERFACE xt_redist_collection_custom_new
    MODULE PROCEDURE xt_redist_collection_new_a_3i_cfg
    MODULE PROCEDURE xt_redist_collection_new_a_2i_cfg
    MODULE PROCEDURE xt_redist_collection_new_a_i_cfg
  END INTERFACE xt_redist_collection_custom_new

  INTERFACE xt_redist_p2p_ext_new
    MODULE PROCEDURE xt_redist_p2p_ext_new_i2_a1d_i2_a1d
    MODULE PROCEDURE xt_redist_p2p_ext_new_i2_a1d_i2_a1d_cfg
    MODULE PROCEDURE xt_redist_p2p_ext_new_i4_a1d_i4_a1d
    MODULE PROCEDURE xt_redist_p2p_ext_new_i4_a1d_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_p2p_ext_new_i8_a1d_i8_a1d
    MODULE PROCEDURE xt_redist_p2p_ext_new_i8_a1d_i8_a1d_cfg
    MODULE PROCEDURE xt_redist_p2p_ext_new_a1d_a1d
    MODULE PROCEDURE xt_redist_p2p_ext_new_a1d_a1d_cfg
  END INTERFACE xt_redist_p2p_ext_new

  INTERFACE xt_redist_p2p_ext_custom_new
    MODULE PROCEDURE xt_redist_p2p_ext_new_i2_a1d_i2_a1d_cfg
    MODULE PROCEDURE xt_redist_p2p_ext_new_i4_a1d_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_p2p_ext_new_i8_a1d_i8_a1d_cfg
    MODULE PROCEDURE xt_redist_p2p_ext_new_a1d_a1d_cfg
  END INTERFACE xt_redist_p2p_ext_custom_new

  INTERFACE xt_redist_repeat_new
    MODULE PROCEDURE xt_redist_repeat_new_i4_a1d
    MODULE PROCEDURE xt_redist_repeat_new_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_repeat_new_a1d
    MODULE PROCEDURE xt_redist_repeat_new_a1d_cfg
    MODULE PROCEDURE xt_redist_repeat_asym_new_i4_a1d
    MODULE PROCEDURE xt_redist_repeat_asym_new_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_repeat_asym_new_a1d
    MODULE PROCEDURE xt_redist_repeat_asym_new_a1d_cfg
  END INTERFACE xt_redist_repeat_new

  INTERFACE xt_redist_repeat_custom_new
    MODULE PROCEDURE xt_redist_repeat_new_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_repeat_new_a1d_cfg
    MODULE PROCEDURE xt_redist_repeat_asym_new_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_repeat_asym_new_a1d_cfg
  END INTERFACE xt_redist_repeat_custom_new

  INTERFACE xt_redist_single_array_base_new
    MODULE PROCEDURE xt_redist_single_array_base_new_i2_a1d_i2_a1d
    MODULE PROCEDURE xt_redist_single_array_base_new_i2_a1d_i2_a1d_cfg
    MODULE PROCEDURE xt_redist_single_array_base_new_i4_a1d_i4_a1d
    MODULE PROCEDURE xt_redist_single_array_base_new_i4_a1d_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_single_array_base_new_i8_a1d_i8_a1d
    MODULE PROCEDURE xt_redist_single_array_base_new_i8_a1d_i8_a1d_cfg
    MODULE PROCEDURE xt_redist_single_array_base_new_a1d_a1d
    MODULE PROCEDURE xt_redist_single_array_base_new_a1d_a1d_cfg
  END INTERFACE xt_redist_single_array_base_new

  INTERFACE xt_redist_single_array_base_custom_new
    MODULE PROCEDURE xt_redist_single_array_base_new_i2_a1d_i2_a1d_cfg
    MODULE PROCEDURE xt_redist_single_array_base_new_i4_a1d_i4_a1d_cfg
    MODULE PROCEDURE xt_redist_single_array_base_new_i8_a1d_i8_a1d_cfg
    MODULE PROCEDURE xt_redist_single_array_base_new_a1d_a1d_cfg
  END INTERFACE xt_redist_single_array_base_custom_new

  PUBLIC :: xt_redist_c2f, xt_redist_f2c, xt_is_null, xt_redist_copy, &
       xt_redist_delete, xt_redist_s_exchange1, xt_redist_s_exchange, &
       xt_redist_collection_static_new, xt_redist_collection_static_custom_new,&
       xt_redist_collection_new, xt_redist_collection_custom_new, &
       xt_redist_repeat_new, xt_redist_repeat_custom_new, &
       xt_redist_get_mpi_comm, xt_redist_p2p_ext_new, &
       xt_redist_a_exchange1, xt_redist_a_exchange, &
       xt_redist_single_array_base_new, xt_redist_single_array_base_custom_new,&
       xt_redist_get_send_mpi_datatype, xt_redist_get_recv_mpi_datatype, &
       xt_redist_get_num_send_msg, xt_redist_get_num_recv_msg
  PUBLIC :: xt_redist_p2p_new, xt_redist_p2p_custom_new
  PUBLIC :: xt_redist_p2p_off_new, xt_redist_p2p_off_custom_new
  PUBLIC :: xt_redist_p2p_blocks_new, xt_redist_p2p_blocks_custom_new
  PUBLIC :: xt_redist_p2p_blocks_off_new, xt_redist_p2p_blocks_off_custom_new
  CHARACTER(len=*), PARAMETER :: filename = 'xt_redist_f.f90'
CONTAINS

  FUNCTION xt_redist_is_null(redist) RESULT(p)
    TYPE(xt_redist), INTENT(in) :: redist
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(redist%cptr)
  END FUNCTION xt_redist_is_null

  FUNCTION xt_redist_c2f(redist) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: redist
    TYPE(xt_redist) :: p
    p%cptr = redist
  END FUNCTION xt_redist_c2f

  FUNCTION xt_redist_copy(redist) RESULT(redist_copy)
    TYPE(xt_redist), INTENT(in) :: redist
    TYPE(xt_redist) :: redist_copy
    INTERFACE
      FUNCTION xt_redist_copy_c(redist) BIND(C, name='xt_redist_copy')
        IMPORT:: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: redist
        TYPE(c_ptr) :: xt_redist_copy_c
      END FUNCTION xt_redist_copy_c
    END INTERFACE
    redist_copy%cptr = xt_redist_copy_c(redist%cptr)
  END FUNCTION xt_redist_copy

  SUBROUTINE xt_redist_delete_1(redist)
    TYPE(xt_redist), INTENT(inout) :: redist
    CALL xt_redist_delete_c(redist%cptr)
    redist%cptr = c_null_ptr
  END SUBROUTINE xt_redist_delete_1

  SUBROUTINE xt_redist_delete_a1d(redists)
    TYPE(xt_redist), INTENT(inout) :: redists(:)
    INTEGER :: i, n
    n = SIZE(redists)
    DO i = 1, n
      CALL xt_redist_delete_c(redists(i)%cptr)
      redists(i)%cptr = c_null_ptr
    END DO
  END SUBROUTINE xt_redist_delete_a1d

  FUNCTION xt_redist_get_num_send_msg(redist) RESULT(num_send_msg)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER :: num_send_msg
    INTEGER(c_int) :: n
    n = xt_redist_get_num_send_msg_c(redist%cptr)
    IF (n > HUGE(num_send_msg) .OR. n < -HUGE(num_send_msg)) &
         CALL xt_abort("num_send_msg out of bounds", filename, __LINE__)
    num_send_msg = INT(n)
  END FUNCTION xt_redist_get_num_send_msg

  FUNCTION xt_redist_get_num_recv_msg(redist) RESULT(num_recv_msg)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER :: num_recv_msg
    INTEGER(c_int) :: n
    n = xt_redist_get_num_recv_msg_c(redist%cptr)
    IF (n > HUGE(num_recv_msg) .OR. n < -HUGE(num_recv_msg)) &
         CALL xt_abort("num_recv_msg out of bounds", filename, __LINE__)
    num_recv_msg = INT(n)
  END FUNCTION xt_redist_get_num_recv_msg

  SUBROUTINE xt_redist_s_exchange1(redist, src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    TYPE(c_ptr), INTENT(in) :: src_data_cptr, dst_data_cptr
    INTERFACE
      SUBROUTINE xt_redist_s_exchange1_c(redist, src_data_cptr, dst_data_cptr) &
           BIND(C, name='xt_redist_s_exchange1')
        IMPORT:: c_ptr
        TYPE(c_ptr), VALUE, INTENT(in) :: redist
        TYPE(c_ptr), VALUE :: src_data_cptr, dst_data_cptr
      END SUBROUTINE xt_redist_s_exchange1_c
    END INTERFACE
    CALL xt_redist_s_exchange1_c(redist%cptr, src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange1

  SUBROUTINE xt_redist_a_exchange1(redist, src_data_cptr, dst_data_cptr, &
       request)
    TYPE(xt_redist), INTENT(in) :: redist
    TYPE(c_ptr) :: src_data_cptr, dst_data_cptr
    TYPE(xt_request), INTENT(out) :: request
    INTERFACE
      SUBROUTINE xt_redist_a_exchange1_c(redist, src_data_cptr, &
           dst_data_cptr, request_c) BIND(C, name='xt_redist_a_exchange1')
        IMPORT:: c_ptr, xt_request
        TYPE(c_ptr), VALUE, INTENT(in) :: redist
        TYPE(c_ptr), VALUE :: src_data_cptr, dst_data_cptr
        TYPE(xt_request), INTENT(out) :: request_c
      END SUBROUTINE xt_redist_a_exchange1_c
    END INTERFACE
    CALL xt_redist_a_exchange1_c(redist%cptr, src_data_cptr, &
         & dst_data_cptr, request)
  END SUBROUTINE xt_redist_a_exchange1

  SUBROUTINE xt_redist_s_exchange_a1d(redist, src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(:), dst_data_cptr(:)
    INTEGER :: n
    INTEGER(c_int) :: num_ptr_c
#if __INTEL_COMPILER == 1910 && __INTEL_COMPILER_UPDATE == 2
    TYPE(c_ptr) :: contig_src_data_cptr(SIZE(src_data_cptr)), &
         contig_dst_data_cptr(SIZE(dst_data_cptr))
    contig_src_data_cptr = src_data_cptr
    contig_dst_data_cptr = dst_data_cptr
#endif
    n = SIZE(src_data_cptr)
    IF (n /= SIZE(dst_data_cptr) .OR. n > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(n, c_int)
#if __INTEL_COMPILER == 1910 && __INTEL_COMPILER_UPDATE == 2
    CALL xt_redist_s_exchange_c(redist%cptr, num_ptr_c, &
         contig_src_data_cptr, contig_dst_data_cptr)
#else
    CALL xt_redist_s_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr)
#endif
  END SUBROUTINE xt_redist_s_exchange_a1d

  SUBROUTINE xt_redist_a_exchange_a1d(redist, src_data_cptr, dst_data_cptr, &
       request)
    TYPE(xt_redist), INTENT(in) :: redist
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(:), dst_data_cptr(:)
    TYPE(xt_request), INTENT(out) :: request
    INTEGER :: num_ptr
    INTEGER(c_int) :: num_ptr_c
#if __INTEL_COMPILER == 1910 && __INTEL_COMPILER_UPDATE == 2
    TYPE(c_ptr) :: contig_src_data_cptr(SIZE(src_data_cptr)), &
         contig_dst_data_cptr(SIZE(dst_data_cptr))
    contig_src_data_cptr = src_data_cptr
    contig_dst_data_cptr = dst_data_cptr
#endif
    num_ptr = SIZE(src_data_cptr)
    IF (num_ptr /= SIZE(dst_data_cptr) .OR. num_ptr > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(num_ptr, c_int)
#if __INTEL_COMPILER == 1910 && __INTEL_COMPILER_UPDATE == 2
    CALL xt_redist_a_exchange_c(redist%cptr, num_ptr_c, &
         contig_src_data_cptr, contig_dst_data_cptr, request)
#else
    CALL xt_redist_a_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr, request)
#endif
  END SUBROUTINE xt_redist_a_exchange_a1d

  SUBROUTINE xt_redist_s_exchange_i2_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i2), INTENT(in) :: num_ptr
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    INTEGER(c_int) :: num_ptr_c
    IF (num_ptr < 0_i2) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(num_ptr, c_int)
    CALL xt_redist_s_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_i2_a1d

  SUBROUTINE xt_redist_a_exchange_i2_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr, request)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i2), INTENT(in) :: num_ptr
    TYPE(xt_request), INTENT(out) :: request
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    INTEGER(c_int) :: num_ptr_c
    IF (num_ptr < 0_i2) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(num_ptr, c_int)
    CALL xt_redist_a_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr, request)
  END SUBROUTINE xt_redist_a_exchange_i2_a1d

  SUBROUTINE xt_redist_s_exchange_i4_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i4), INTENT(in) :: num_ptr
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    INTEGER(c_int) :: num_ptr_c
    IF (num_ptr < 0_i4 .OR. num_ptr > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(num_ptr, c_int)
    CALL xt_redist_s_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_i4_a1d

  SUBROUTINE xt_redist_a_exchange_i4_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr, request)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i4), INTENT(in) :: num_ptr
    TYPE(xt_request), INTENT(out) :: request
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    INTEGER(c_int) :: num_ptr_c
    IF (num_ptr < 0_i4 .OR. num_ptr > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(num_ptr, c_int)
    CALL xt_redist_a_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr, request)
  END SUBROUTINE  xt_redist_a_exchange_i4_a1d

  SUBROUTINE xt_redist_s_exchange_i8_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i8), INTENT(in) :: num_ptr
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    INTEGER(c_int) :: num_ptr_c
    IF (num_ptr < 0_i8 .OR. num_ptr > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(num_ptr, c_int)
    CALL xt_redist_s_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_i8_a1d

  SUBROUTINE xt_redist_a_exchange_i8_a1d(redist, num_ptr, &
       src_data_cptr, dst_data_cptr, request)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(i8), INTENT(in) :: num_ptr
    TYPE(xt_request), INTENT(out) :: request
    TYPE(c_ptr), INTENT(in) :: src_data_cptr(num_ptr), dst_data_cptr(num_ptr)
    INTEGER(c_int) :: num_ptr_c
    IF (num_ptr < 0_i8 .OR. num_ptr > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of pointers", filename, __LINE__)
    num_ptr_c = INT(num_ptr, c_int)
    CALL xt_redist_a_exchange_c(redist%cptr, num_ptr_c, &
         src_data_cptr, dst_data_cptr, request)
  END SUBROUTINE xt_redist_a_exchange_i8_a1d

  FUNCTION xt_redist_p2p_new(xmap, datatype) RESULT(res)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_new_f(xmap, datatype) &
           BIND(C, name='xt_redist_p2p_new_f') RESULT(res_ptr)
        IMPORT:: xt_xmap, xt_redist, c_int, c_ptr, xt_mpi_fint_kind
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_new_f
    END INTERFACE
    res%cptr = xt_redist_p2p_new_f(xmap, datatype)
  END FUNCTION xt_redist_p2p_new

  FUNCTION xt_redist_p2p_custom_new(xmap, datatype, config) RESULT(res)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_custom_new_f(xmap, datatype, config) &
           BIND(C, name='xt_redist_p2p_custom_new_f') RESULT(res_ptr)
        IMPORT:: xt_xmap, xt_redist, c_ptr, xt_mpi_fint_kind, xt_config
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(xt_config), INTENT(in) :: config
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_custom_new_f
    END INTERFACE
    res%cptr = xt_redist_p2p_custom_new_f(xmap, datatype, config)
  END FUNCTION xt_redist_p2p_custom_new

  FUNCTION xt_redist_p2p_off_new(xmap, src_offsets, dst_offsets, datatype) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: src_offsets(*)
    INTEGER, INTENT(in) :: dst_offsets(*)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_off_new_f(xmap, src_offsets, dst_offsets, &
           datatype) BIND(C, name='xt_redist_p2p_off_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, c_ptr, xt_mpi_fint_kind
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: src_offsets(*), dst_offsets(*)
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_off_new_f
    END INTERFACE
    res%cptr = xt_redist_p2p_off_new_f(xmap, src_offsets, dst_offsets, datatype)
  END FUNCTION xt_redist_p2p_off_new

  FUNCTION xt_redist_p2p_off_custom_new(xmap, src_offsets, dst_offsets, &
       datatype, config) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER, INTENT(in) :: src_offsets(*), dst_offsets(*)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_off_custom_new_f(xmap, src_offsets, dst_offsets, &
           datatype, config) BIND(C, name='xt_redist_p2p_off_custom_new_f') &
           RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, c_ptr, xt_mpi_fint_kind, xt_config
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(xt_mpi_fint_kind), INTENT(in) :: src_offsets(*), dst_offsets(*)
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(xt_config), INTENT(in) :: config
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_off_custom_new_f
    END INTERFACE
    res%cptr = xt_redist_p2p_off_custom_new_f(xmap, src_offsets, dst_offsets, &
         datatype, config)
  END FUNCTION xt_redist_p2p_off_custom_new

  FUNCTION xt_redist_p2p_blocks_new(xmap, src_block_sizes, src_block_num, &
       &                                  dst_block_sizes, dst_block_num, &
       &                                  datatype) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_sizes(*), src_block_num, &
         dst_block_sizes(*), dst_block_num
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_new_f(xmap, &
           &                              src_block_sizes, src_block_num, &
           &                              dst_block_sizes, dst_block_num, &
           &                              datatype) &
           BIND(C, name='xt_redist_p2p_blocks_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_mpi_fint_kind, xt_redist, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num, dst_block_num
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*), dst_block_sizes(*)
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_new_f
    END INTERFACE
    res%cptr = xt_redist_p2p_blocks_new_f(xmap, &
         &                                src_block_sizes, src_block_num, &
         &                                dst_block_sizes, dst_block_num, &
         &                                datatype)
  END FUNCTION xt_redist_p2p_blocks_new

  FUNCTION xt_redist_p2p_blocks_custom_new(xmap, src_block_sizes, &
       &                                   src_block_num, &
       &                                   dst_block_sizes, dst_block_num, &
       &                                   datatype, config) &
       RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_sizes(*), src_block_num, &
         dst_block_sizes(*), dst_block_num
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_custom_new_f(xmap, &
           &                              src_block_sizes, src_block_num, &
           &                              dst_block_sizes, dst_block_num, &
           &                              datatype, config) &
           BIND(C, name='xt_redist_p2p_blocks_custom_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_mpi_fint_kind, xt_redist, c_int, c_ptr, xt_config
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num, dst_block_num
        INTEGER(c_int), INTENT(in) :: src_block_sizes(*), dst_block_sizes(*)
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(xt_config), INTENT(in) :: config
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_custom_new_f
    END INTERFACE
    res%cptr = xt_redist_p2p_blocks_custom_new_f(xmap, &
         &  src_block_sizes, src_block_num, &
         &  dst_block_sizes, dst_block_num, datatype, config)
  END FUNCTION xt_redist_p2p_blocks_custom_new

  FUNCTION xt_redist_p2p_blocks_off_new(xmap, src_block_offsets, &
       src_block_sizes, src_block_num, &
       dst_block_offsets, dst_block_sizes, dst_block_num, &
       datatype) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_offsets(*)
    INTEGER(c_int), INTENT(in) :: src_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
    INTEGER(c_int), INTENT(in) :: dst_block_offsets(*)
    INTEGER(c_int), INTENT(in) :: dst_block_sizes(*)
    INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_off_new_f(xmap, src_block_offsets, &
           src_block_sizes, src_block_num, &
           dst_block_offsets, dst_block_sizes, dst_block_num, &
           datatype) BIND(C, name='xt_redist_p2p_blocks_off_new_f') &
           RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, xt_mpi_fint_kind, c_int, c_ptr
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_offsets(*), src_block_sizes(*),&
             dst_block_offsets(*), dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num, dst_block_num
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_off_new_f
    END INTERFACE
    res%cptr = xt_redist_p2p_blocks_off_new_f(xmap, src_block_offsets, &
         src_block_sizes, src_block_num, &
         dst_block_offsets, dst_block_sizes, dst_block_num, datatype)
  END FUNCTION xt_redist_p2p_blocks_off_new

  FUNCTION xt_redist_p2p_blocks_off_custom_new(xmap, src_block_offsets, &
       src_block_sizes, src_block_num, &
       dst_block_offsets, dst_block_sizes, dst_block_num, &
       datatype, config) RESULT(res)
    IMPLICIT NONE
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(c_int), INTENT(in) :: src_block_offsets(*), src_block_sizes(*), &
         dst_block_offsets(*), dst_block_sizes(*)
    INTEGER(c_int), INTENT(in) :: src_block_num, dst_block_num
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTERFACE
      FUNCTION xt_redist_p2p_blocks_off_custom_new_f(xmap, src_block_offsets, &
           src_block_sizes, src_block_num, &
           dst_block_offsets, dst_block_sizes, dst_block_num, &
           datatype, config) &
           BIND(C, name='xt_redist_p2p_blocks_off_custom_new_f') RESULT(res_ptr)
        IMPORT :: xt_xmap, xt_redist, xt_mpi_fint_kind, c_int, c_ptr, xt_config
        IMPLICIT NONE
        TYPE(xt_xmap), INTENT(in) :: xmap
        INTEGER(c_int), INTENT(in) :: src_block_offsets(*), src_block_sizes(*),&
             dst_block_offsets(*), dst_block_sizes(*)
        INTEGER(c_int), VALUE, INTENT(in) :: src_block_num
        INTEGER(c_int), VALUE, INTENT(in) :: dst_block_num
        INTEGER(xt_mpi_fint_kind), VALUE, INTENT(in) :: datatype
        TYPE(xt_config), INTENT(in) :: config
        TYPE(c_ptr) :: res_ptr
      END FUNCTION xt_redist_p2p_blocks_off_custom_new_f
    END INTERFACE

    res%cptr = xt_redist_p2p_blocks_off_custom_new_f(xmap, &
         src_block_offsets, src_block_sizes, src_block_num, &
         dst_block_offsets, dst_block_sizes, dst_block_num, &
         datatype, config)
  END FUNCTION xt_redist_p2p_blocks_off_custom_new

  FUNCTION xt_redist_collection_static_new_a_i_2ak_i(redists, num_redists, &
       src_displacements, dst_displacements, comm) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(*)
    INTEGER, INTENT(in) :: num_redists, comm
    INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(*), &
         dst_displacements(*)
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_redists_c
    num_redists_c = INT(num_redists, c_int)
    res%cptr = xt_redist_collection_static_new_f(redists, &
         num_redists_c, src_displacements, dst_displacements, comm)
  END FUNCTION xt_redist_collection_static_new_a_i_2ak_i

  FUNCTION xt_redist_collection_static_new_a_i_2ak_i_cfg(redists, num_redists, &
       src_displacements, dst_displacements, comm, config) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(*)
    INTEGER, INTENT(in) :: num_redists, comm
    INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(*), &
         dst_displacements(*)
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_redists_c
    num_redists_c = INT(num_redists, c_int)
    res%cptr = xt_redist_collection_static_custom_new_f(redists, &
         num_redists_c, src_displacements, dst_displacements, comm, config)
  END FUNCTION xt_redist_collection_static_new_a_i_2ak_i_cfg

  FUNCTION xt_redist_collection_static_new_a_2ak_i(redists, &
       src_displacements, dst_displacements, comm) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(:)
    INTEGER, INTENT(in) :: comm
    INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(:), &
         dst_displacements(:)
    TYPE(xt_redist) :: res
    INTEGER :: num_redists
    num_redists = SIZE(redists)
    IF (num_redists /= SIZE(src_displacements) &
         .OR. num_redists /= SIZE(dst_displacements)) &
         CALL xt_abort("invalid number of redists", filename, __LINE__)
    res%cptr = xt_redist_collection_static_new_f(redists, &
         num_redists, src_displacements, dst_displacements, comm)
  END FUNCTION xt_redist_collection_static_new_a_2ak_i

  FUNCTION xt_redist_collection_static_new_a_2ak_i_cfg(redists, &
       src_displacements, dst_displacements, comm, config) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(:)
    INTEGER, INTENT(in) :: comm
    INTEGER(mpi_address_kind), INTENT(in) :: src_displacements(:), &
         dst_displacements(:)
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTEGER :: num_redists
    num_redists = SIZE(redists)
    IF (num_redists /= SIZE(src_displacements) &
         .OR. num_redists /= SIZE(dst_displacements)) &
         CALL xt_abort("invalid number of redists", filename, __LINE__)
    res%cptr = xt_redist_collection_static_custom_new_f(redists, &
         num_redists, src_displacements, dst_displacements, comm, config)
  END FUNCTION xt_redist_collection_static_new_a_2ak_i_cfg

  FUNCTION xt_redist_collection_new_a_3i(redists, num_redists, cache_size, &
       comm) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(*)
    INTEGER, INTENT(in) :: num_redists, cache_size, comm
    TYPE(xt_redist) :: res
    res%cptr = xt_redist_collection_new_f(redists, num_redists, &
         cache_size, comm)
  END FUNCTION xt_redist_collection_new_a_3i

  FUNCTION xt_redist_collection_new_a_3i_cfg(redists, num_redists, cache_size, &
       comm, config) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(*)
    INTEGER, INTENT(in) :: num_redists, cache_size, comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    res%cptr = xt_redist_collection_custom_new_f(redists, num_redists,&
         cache_size, comm, config)
  END FUNCTION xt_redist_collection_new_a_3i_cfg

  FUNCTION xt_redist_collection_new_a_2i(redists, cache_size, comm) &
       RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(:)
    INTEGER, INTENT(in) :: cache_size, comm
    TYPE(xt_redist) :: res
    INTEGER :: num_redists
    num_redists = SIZE(redists)
    res%cptr = xt_redist_collection_new_f(redists, &
         num_redists, cache_size, comm)
  END FUNCTION xt_redist_collection_new_a_2i

  FUNCTION xt_redist_collection_new_a_2i_cfg(redists, cache_size, comm, &
       config) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(:)
    INTEGER, INTENT(in) :: cache_size, comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTEGER :: num_redists
    num_redists = SIZE(redists)
    res%cptr = xt_redist_collection_custom_new_f(redists, &
         num_redists, cache_size, comm, config)
  END FUNCTION xt_redist_collection_new_a_2i_cfg

  FUNCTION xt_redist_collection_new_a_i(redists, comm) &
       RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(:)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_redist) :: res

    res = xt_redist_collection_new_a_3i(redists, SIZE(redists), -1, comm)
  END FUNCTION xt_redist_collection_new_a_i

  FUNCTION xt_redist_collection_new_a_i_cfg(redists, comm, config) &
       RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redists(:)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res

    res = xt_redist_collection_new_a_3i_cfg(redists, SIZE(redists), &
         &-1, comm, config)
  END FUNCTION xt_redist_collection_new_a_i_cfg


  FUNCTION xt_redist_repeat_new_i4_a1d(redist, src_extent, dst_extent, &
       num_repetitions, displacements) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(i4), INTENT(in) :: num_repetitions
    INTEGER(c_int), INTENT(in) :: displacements(num_repetitions)
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_repetitions_c
    num_repetitions_c = INT(num_repetitions, c_int)
    res%cptr = xt_redist_repeat_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions_c, displacements)
  END FUNCTION xt_redist_repeat_new_i4_a1d

  FUNCTION xt_redist_repeat_new_i4_a1d_cfg(redist, src_extent, dst_extent, &
       num_repetitions, displacements, config) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(i4), INTENT(in) :: num_repetitions
    INTEGER(c_int), INTENT(in) :: displacements(num_repetitions)
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_repetitions_c
    num_repetitions_c = INT(num_repetitions, c_int)
    res%cptr = xt_redist_repeat_custom_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions_c, displacements, &
         xt_config_f2c(config))
  END FUNCTION xt_redist_repeat_new_i4_a1d_cfg

  FUNCTION xt_redist_repeat_new_a1d(redist, src_extent, dst_extent, &
       displacements) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(c_int), INTENT(in) :: displacements(:)
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_repetitions_c
    num_repetitions_c = INT(SIZE(displacements), c_int)
    res%cptr = xt_redist_repeat_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions_c, displacements)
  END FUNCTION xt_redist_repeat_new_a1d

  FUNCTION xt_redist_repeat_new_a1d_cfg(redist, src_extent, dst_extent, &
       displacements, config) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(c_int), INTENT(in) :: displacements(:)
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_repetitions_c
    num_repetitions_c = INT(SIZE(displacements), c_int)
    res%cptr = xt_redist_repeat_custom_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions_c, displacements, &
         xt_config_f2c(config))
  END FUNCTION xt_redist_repeat_new_a1d_cfg

  FUNCTION xt_redist_repeat_asym_new_i4_a1d(redist, src_extent, dst_extent, &
       num_repetitions, src_displacements, dst_displacements) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(i4), INTENT(in) :: num_repetitions
    INTEGER(c_int), INTENT(in) :: src_displacements(num_repetitions), &
         & dst_displacements(num_repetitions)
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_repetitions_c
    num_repetitions_c = INT(num_repetitions, c_int)
    res%cptr = xt_redist_repeat_asym_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions_c, src_displacements, &
         dst_displacements)
  END FUNCTION xt_redist_repeat_asym_new_i4_a1d

  FUNCTION xt_redist_repeat_asym_new_i4_a1d_cfg(redist, src_extent, dst_extent,&
       num_repetitions, src_displacements, dst_displacements, config) &
       RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(i4), INTENT(in) :: num_repetitions
    INTEGER(c_int), INTENT(in) :: src_displacements(num_repetitions), &
         & dst_displacements(num_repetitions)
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTEGER(c_int) :: num_repetitions_c
    num_repetitions_c = INT(num_repetitions, c_int)
    res%cptr = xt_redist_repeat_asym_custom_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions_c, src_displacements, &
         dst_displacements, xt_config_f2c(config))
  END FUNCTION xt_redist_repeat_asym_new_i4_a1d_cfg

  FUNCTION xt_redist_repeat_asym_new_a1d(redist, src_extent, dst_extent, &
       src_displacements, dst_displacements) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(c_int), INTENT(in) :: src_displacements(:), dst_displacements(:)
    TYPE(xt_redist) :: res
    INTEGER(i4) :: num_repetitions
    num_repetitions = SIZE(src_displacements)
    IF (num_repetitions /= SIZE(dst_displacements)) &
         CALL xt_abort("inequal size for src and dst displacements", &
         filename, __LINE__)
    res%cptr = xt_redist_repeat_asym_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions, src_displacements, &
         dst_displacements)
  END FUNCTION xt_redist_repeat_asym_new_a1d

  FUNCTION xt_redist_repeat_asym_new_a1d_cfg(redist, src_extent, dst_extent, &
       src_displacements, dst_displacements, config) RESULT(res)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER(mpi_address_kind), INTENT(in) :: src_extent, dst_extent
    INTEGER(c_int), INTENT(in) :: src_displacements(:), dst_displacements(:)
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: res
    INTEGER(i4) :: num_repetitions
    num_repetitions = SIZE(src_displacements)
    IF (num_repetitions /= SIZE(dst_displacements)) &
         CALL xt_abort("inequal size for src and dst displacements", &
         filename, __LINE__)
    res%cptr = xt_redist_repeat_asym_custom_new_c(redist%cptr, &
         src_extent, dst_extent, num_repetitions, src_displacements, &
         dst_displacements, xt_config_f2c(config))
  END FUNCTION xt_redist_repeat_asym_new_a1d_cfg

  FUNCTION xt_redist_p2p_ext_new_i2_a1d_i2_a1d(xmap, num_src_ext, src_extents, &
       num_dst_ext, dst_extents, datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i2), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    IF (num_src_ext < 0_i2 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i2 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_new_c2f(xmap, &
         num_src_ext_c, src_extents, num_dst_ext_c, dst_extents, datatype)
  END FUNCTION xt_redist_p2p_ext_new_i2_a1d_i2_a1d

  FUNCTION xt_redist_p2p_ext_new_i2_a1d_i2_a1d_cfg(xmap, num_src_ext, &
       src_extents, num_dst_ext, dst_extents, datatype, config) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i2), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    IF (num_src_ext < 0_i2 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i2 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_custom_new_c2f(xmap, num_src_ext_c, &
         src_extents, num_dst_ext_c, dst_extents, datatype, config)
  END FUNCTION xt_redist_p2p_ext_new_i2_a1d_i2_a1d_cfg

  FUNCTION xt_redist_p2p_ext_new_i4_a1d_i4_a1d(xmap, num_src_ext, src_extents, &
       num_dst_ext, dst_extents, datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i4), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    IF (num_src_ext < 0_i4 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i4 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_new_c2f(xmap, num_src_ext_c, &
         src_extents, num_dst_ext_c, dst_extents, datatype)
  END FUNCTION xt_redist_p2p_ext_new_i4_a1d_i4_a1d

  FUNCTION xt_redist_p2p_ext_new_i4_a1d_i4_a1d_cfg(xmap, num_src_ext, &
       src_extents, num_dst_ext, dst_extents, datatype, config) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i4), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    IF (num_src_ext < 0_i4 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i4 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_custom_new_c2f(xmap, &
         num_src_ext_c, src_extents, num_dst_ext_c, dst_extents, datatype, &
         config)
  END FUNCTION xt_redist_p2p_ext_new_i4_a1d_i4_a1d_cfg

  FUNCTION xt_redist_p2p_ext_new_i8_a1d_i8_a1d(xmap, num_src_ext, src_extents, &
       num_dst_ext, dst_extents, datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i8), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    IF (num_src_ext < 0_i8 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i8 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_new_c2f(xmap, num_src_ext_c, &
         src_extents, num_dst_ext_c, dst_extents, datatype)
  END FUNCTION xt_redist_p2p_ext_new_i8_a1d_i8_a1d

  FUNCTION xt_redist_p2p_ext_new_i8_a1d_i8_a1d_cfg(xmap, num_src_ext, &
       src_extents, num_dst_ext, dst_extents, datatype, config) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    INTEGER(i8), INTENT(in) :: num_src_ext, num_dst_ext
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(num_src_ext), &
         dst_extents(num_dst_ext)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    IF (num_src_ext < 0_i8 .OR. num_src_ext > HUGE(1_c_int) &
         .OR. num_dst_ext < 0_i8 .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_custom_new_c2f(xmap, num_src_ext_c, &
         src_extents, num_dst_ext_c, dst_extents, datatype, config)
  END FUNCTION xt_redist_p2p_ext_new_i8_a1d_i8_a1d_cfg

  FUNCTION xt_redist_p2p_ext_new_a1d_a1d(xmap, src_extents, dst_extents, &
       datatype) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(:), dst_extents(:)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_redist) :: redist
    INTEGER :: num_src_ext, num_dst_ext
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    num_src_ext = SIZE(src_extents)
    num_dst_ext = SIZE(dst_extents)
    IF (num_src_ext > HUGE(1_c_int) .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_new_c2f(xmap, num_src_ext_c, &
         src_extents, num_dst_ext_c, dst_extents, datatype)
  END FUNCTION xt_redist_p2p_ext_new_a1d_a1d

  FUNCTION xt_redist_p2p_ext_new_a1d_a1d_cfg(xmap, src_extents, dst_extents, &
       datatype, config) RESULT(redist)
    TYPE(xt_xmap), INTENT(in) :: xmap
    TYPE(xt_offset_ext), INTENT(in) :: src_extents(:), dst_extents(:)
    INTEGER, INTENT(in) :: datatype
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist
    INTEGER :: num_src_ext, num_dst_ext
    INTEGER(c_int) :: num_src_ext_c, num_dst_ext_c
    num_src_ext = SIZE(src_extents)
    num_dst_ext = SIZE(dst_extents)
    IF (num_src_ext > HUGE(1_c_int) .OR. num_dst_ext > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of extents", filename, __LINE__)
    num_src_ext_c = INT(num_src_ext, c_int)
    num_dst_ext_c = INT(num_dst_ext, c_int)
    redist%cptr = xt_redist_p2p_ext_custom_new_c2f(xmap, num_src_ext_c, &
         src_extents, num_dst_ext_c, dst_extents, datatype, config)
  END FUNCTION xt_redist_p2p_ext_new_a1d_a1d_cfg

  FUNCTION xt_redist_single_array_base_new_i2_a1d_i2_a1d( &
       nsend, nrecv, send_msgs, recv_msgs, comm) RESULT(redist)
    INTEGER(i2), INTENT(in) :: nsend, nrecv
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: &
         send_msgs(nsend), recv_msgs(nrecv)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_redist) :: redist

    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p

    IF (nsend < 0_i2 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0_i2 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    IF (nsend > 0_i2) THEN
      send_msgs_p = C_LOC(send_msgs)
    ELSE
      send_msgs_p = c_null_ptr
    END IF
    IF (nrecv > 0_i2) THEN
      recv_msgs_p = C_LOC(recv_msgs)
    ELSE
      recv_msgs_p = c_null_ptr
    END IF
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist%cptr = xt_redist_single_array_base_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm)
  END FUNCTION xt_redist_single_array_base_new_i2_a1d_i2_a1d

  FUNCTION xt_redist_single_array_base_new_i2_a1d_i2_a1d_cfg( &
       nsend, nrecv, send_msgs, recv_msgs, comm, config) RESULT(redist)
    INTEGER(i2), INTENT(in) :: nsend, nrecv
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: &
         send_msgs(nsend), recv_msgs(nrecv)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist

    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p

    IF (nsend < 0_i2 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0_i2 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    IF (nsend > 0_i2) THEN
      send_msgs_p = C_LOC(send_msgs)
    ELSE
      send_msgs_p = c_null_ptr
    END IF
    IF (nrecv > 0_i2) THEN
      recv_msgs_p = C_LOC(recv_msgs)
    ELSE
      recv_msgs_p = c_null_ptr
    END IF
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist%cptr = xt_redist_single_array_base_custom_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm, config)
  END FUNCTION xt_redist_single_array_base_new_i2_a1d_i2_a1d_cfg

  FUNCTION xt_redist_single_array_base_new_i4_a1d_i4_a1d( &
       nsend, nrecv, send_msgs, recv_msgs, comm) RESULT(redist)
    INTEGER(i4), INTENT(in) :: nsend, nrecv
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: &
         send_msgs(nsend), recv_msgs(nrecv)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_redist) :: redist

    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p

    IF (nsend < 0_i4 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0_i4 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    IF (nsend > 0_i4) THEN
      send_msgs_p = C_LOC(send_msgs)
    ELSE
      send_msgs_p = c_null_ptr
    END IF
    IF (nrecv > 0_i4) THEN
      recv_msgs_p = C_LOC(recv_msgs)
    ELSE
      recv_msgs_p = c_null_ptr
    END IF
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist = xt_redist_c2f(xt_redist_single_array_base_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm))
  END FUNCTION xt_redist_single_array_base_new_i4_a1d_i4_a1d

  FUNCTION xt_redist_single_array_base_new_i4_a1d_i4_a1d_cfg( &
       nsend, nrecv, send_msgs, recv_msgs, comm, config) RESULT(redist)
    INTEGER(i4), INTENT(in) :: nsend, nrecv
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: &
         send_msgs(nsend), recv_msgs(nrecv)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist

    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p

    IF (nsend < 0_i4 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0_i4 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    IF (nsend > 0_i4) THEN
      send_msgs_p = C_LOC(send_msgs)
    ELSE
      send_msgs_p = c_null_ptr
    END IF
    IF (nrecv > 0_i4) THEN
      recv_msgs_p = C_LOC(recv_msgs)
    ELSE
      recv_msgs_p = c_null_ptr
    END IF
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist%cptr = xt_redist_single_array_base_custom_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm, config)
  END FUNCTION xt_redist_single_array_base_new_i4_a1d_i4_a1d_cfg

  FUNCTION xt_redist_single_array_base_new_i8_a1d_i8_a1d( &
       nsend, nrecv, send_msgs, recv_msgs, comm) RESULT(redist)
    INTEGER(i8), INTENT(in) :: nsend, nrecv
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: &
         send_msgs(nsend), recv_msgs(nrecv)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_redist) :: redist

    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p

    IF (nsend < 0_i8 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0_i8 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    IF (nsend > 0_i8) THEN
      send_msgs_p = C_LOC(send_msgs)
    ELSE
      send_msgs_p = c_null_ptr
    END IF
    IF (nrecv > 0_i8) THEN
      recv_msgs_p = C_LOC(recv_msgs)
    ELSE
      recv_msgs_p = c_null_ptr
    END IF
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist%cptr = xt_redist_single_array_base_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm)
  END FUNCTION xt_redist_single_array_base_new_i8_a1d_i8_a1d

  FUNCTION xt_redist_single_array_base_new_i8_a1d_i8_a1d_cfg( &
       nsend, nrecv, send_msgs, recv_msgs, comm, config) RESULT(redist)
    INTEGER(i8), INTENT(in) :: nsend, nrecv
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: &
         send_msgs(nsend), recv_msgs(nrecv)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist

    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p

    IF (nsend < 0_i8 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0_i8 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    IF (nsend > 0_i8) THEN
      send_msgs_p = C_LOC(send_msgs)
    ELSE
      send_msgs_p = c_null_ptr
    END IF
    IF (nrecv > 0_i8) THEN
      recv_msgs_p = C_LOC(recv_msgs)
    ELSE
      recv_msgs_p = c_null_ptr
    END IF
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist%cptr = xt_redist_single_array_base_custom_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm, config)
  END FUNCTION xt_redist_single_array_base_new_i8_a1d_i8_a1d_cfg

  FUNCTION xt_redist_single_array_base_new_a1d_a1d(send_msgs, recv_msgs, comm) &
       RESULT(redist)
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: send_msgs(:), recv_msgs(:)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_redist) :: redist

    INTEGER :: nsend, nrecv
    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p
    TYPE(xt_redist_msg), ALLOCATABLE, TARGET :: send_msgs_a(:), recv_msgs_a(:)

    nsend = SIZE(send_msgs)
    nrecv = SIZE(recv_msgs)

    CALL msgs_p_arg(send_msgs, send_msgs_a, send_msgs_p)
    CALL msgs_p_arg(recv_msgs, recv_msgs_a, recv_msgs_p)

    IF (nsend < 0 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist%cptr = xt_redist_single_array_base_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm)
  END FUNCTION xt_redist_single_array_base_new_a1d_a1d

  FUNCTION xt_redist_single_array_base_new_a1d_a1d_cfg(send_msgs, recv_msgs, &
       comm, config) RESULT(redist)
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: send_msgs(:), recv_msgs(:)
    INTEGER, INTENT(in) :: comm
    TYPE(xt_config), INTENT(in) :: config
    TYPE(xt_redist) :: redist

    INTEGER :: nsend, nrecv
    INTEGER(c_int) :: nsend_c, nrecv_c
    TYPE(c_ptr) :: send_msgs_p, recv_msgs_p
    TYPE(xt_redist_msg), ALLOCATABLE, TARGET :: send_msgs_a(:), recv_msgs_a(:)

    nsend = SIZE(send_msgs)
    nrecv = SIZE(recv_msgs)

    CALL msgs_p_arg(send_msgs, send_msgs_a, send_msgs_p)
    CALL msgs_p_arg(recv_msgs, recv_msgs_a, recv_msgs_p)

    IF (nsend < 0 .OR. nsend > HUGE(1_c_int) &
         .OR. nrecv < 0 .OR. nrecv > HUGE(1_c_int)) &
         CALL xt_abort("invalid number of send messages", filename, __LINE__)
    nsend_c = INT(nsend, c_int)
    nrecv_c = INT(nrecv, c_int)
    redist%cptr = xt_redist_single_array_base_custom_new_c2f(&
         nsend_c, nrecv_c, send_msgs_p, recv_msgs_p, comm, config)
  END FUNCTION xt_redist_single_array_base_new_a1d_a1d_cfg

  SUBROUTINE msgs_p_arg(msgs, msgs_a, msgs_p)
    TYPE(xt_redist_msg), TARGET, INTENT(in) :: msgs(:)
    TYPE(xt_redist_msg), TARGET, ALLOCATABLE, INTENT(inout) :: msgs_a(:)
    TYPE(c_ptr), INTENT(out) :: msgs_p

    INTEGER :: msgs_size
    LOGICAL :: msgs_is_contiguous
#ifndef HAVE_FC_IS_CONTIGUOUS
    INTERFACE
      FUNCTION xt_redist_msg_contiguous(msgs_a, msgs_b) RESULT(p) &
            BIND(c, name='xt_redist_msg_contiguous')
        IMPORT :: c_int, xt_redist_msg
        TYPE(xt_redist_msg), INTENT(in) :: msgs_a, msgs_b
        INTEGER(c_int) :: p
      END FUNCTION xt_redist_msg_contiguous
    END INTERFACE
#endif

    msgs_size = SIZE(msgs)
    IF (msgs_size > HUGE(1_c_int)) &
         CALL xt_abort('invalid size', filename, __LINE__)
    IF (msgs_size > 0) THEN
      IF (msgs_size > 1) THEN
#ifdef HAVE_FC_IS_CONTIGUOUS
        msgs_is_contiguous = IS_CONTIGUOUS(msgs)
#else
        msgs_is_contiguous = xt_redist_msg_contiguous(msgs(1), msgs(2)) /= 0
#endif
        IF (msgs_is_contiguous) THEN
          XT_SLICE_C_LOC(msgs(1), msgs_p)
        ELSE
          ALLOCATE(msgs_a(msgs_size))
          msgs_a = msgs
          msgs_p = C_LOC(msgs_a)
        END IF
      ELSE
        XT_SLICE_C_LOC(msgs(1), msgs_p)
      END IF
    ELSE
      msgs_p = c_null_ptr
    END IF
  END SUBROUTINE msgs_p_arg

END MODULE xt_redist_base

MODULE xt_redist_rename
  USE xt_redist_base, ONLY: xt_redist_p2p_orig_new => xt_redist_p2p_new, &
       xt_redist_p2p_custom_new
  IMPLICIT NONE
  PRIVATE
  INTERFACE xt_redist_p2p_new
    MODULE PROCEDURE xt_redist_p2p_orig_new
    MODULE PROCEDURE xt_redist_p2p_custom_new
  END INTERFACE xt_redist_p2p_new
  PUBLIC :: xt_redist_p2p_new
END MODULE xt_redist_rename
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
