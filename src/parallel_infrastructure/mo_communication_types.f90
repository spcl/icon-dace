! This module provides communication pattern types.
! for parallel runs
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

#include "crayftn_ptr_fail.inc"
MODULE mo_communication_types
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_kind, ONLY: sp, dp
USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup
USE mo_fortran_tools, ONLY: t_ptr_3d, t_ptr_3d_sp

IMPLICIT NONE

PRIVATE

PUBLIC :: t_comm_pattern, t_comm_pattern_collection, t_p_comm_pattern

  TYPE xfer_list
    INTEGER :: rank
    INTEGER, ALLOCATABLE :: glob_idx(:)
  END TYPE xfer_list
  PUBLIC :: xfer_list

!--------------------------------------------------------------------------------------------------
!
TYPE, ABSTRACT :: t_comm_pattern

  CONTAINS

    PROCEDURE(interface_setup_comm_pattern), DEFERRED :: setup
    PROCEDURE(interface_setup_comm_pattern2), DEFERRED :: setup2
    PROCEDURE(interface_delete_comm_pattern), DEFERRED :: delete
    PROCEDURE(interface_exchange_data_r3d), DEFERRED :: exchange_data_r3d
    PROCEDURE(interface_exchange_data_s3d), DEFERRED :: exchange_data_s3d
    PROCEDURE(interface_exchange_data_i3d), DEFERRED :: exchange_data_i3d
    PROCEDURE(interface_exchange_data_r2d), DEFERRED :: exchange_data_r2d
    PROCEDURE(interface_exchange_data_s2d), DEFERRED :: exchange_data_s2d
    PROCEDURE(interface_exchange_data_i2d), DEFERRED :: exchange_data_i2d
    PROCEDURE(interface_exchange_data_l2d), DEFERRED :: exchange_data_l2d
    PROCEDURE(interface_exchange_data_l3d), DEFERRED :: exchange_data_l3d
    PROCEDURE(interface_exchange_data_mult), DEFERRED :: exchange_data_mult
    PROCEDURE(interface_exchange_data_mult_mixprec), DEFERRED :: exchange_data_mult_mixprec
    PROCEDURE(interface_exchange_data_4de1), DEFERRED :: exchange_data_4de1
    PROCEDURE(interface_get_np_recv), DEFERRED :: get_np_recv
    PROCEDURE(interface_get_np_send), DEFERRED :: get_np_send
    PROCEDURE(interface_get_pelist_recv), DEFERRED :: get_pelist_recv
END TYPE t_comm_pattern

TYPE, ABSTRACT :: t_comm_pattern_collection
  CONTAINS

    PROCEDURE(interface_setup_comm_pattern_collection), DEFERRED :: setup
    PROCEDURE(interface_delete_comm_pattern_collection), DEFERRED :: delete
    PROCEDURE(interface_exchange_data_grf), DEFERRED :: exchange_data_grf
END TYPE t_comm_pattern_collection

TYPE t_p_comm_pattern
  CLASS(t_comm_pattern), POINTER :: p
END TYPE t_p_comm_pattern

ABSTRACT INTERFACE

  SUBROUTINE interface_setup_comm_pattern_collection(pattern_collection, &
                                                     patterns)
    IMPORT t_comm_pattern_collection, t_p_comm_pattern
    CLASS(t_comm_pattern_collection), INTENT(OUT) :: pattern_collection
    TYPE(t_p_comm_pattern), INTENT(IN) :: patterns(:)
  END SUBROUTINE interface_setup_comm_pattern_collection

  SUBROUTINE interface_setup_comm_pattern( &
    p_pat, dst_n_points, dst_owner, dst_global_index, send_glb2loc_index, &
    src_n_points, src_owner, src_global_index, inplace, comm)
    IMPORT t_comm_pattern, t_glb2loc_index_lookup
    CLASS(t_comm_pattern), TARGET, INTENT(OUT) :: p_pat
    INTEGER, INTENT(IN)           :: dst_n_points
    INTEGER, INTENT(IN)           :: dst_owner(:)
    INTEGER, INTENT(IN)           :: dst_global_index(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: send_glb2loc_index
    INTEGER, INTENT(IN)           :: src_n_points
    INTEGER, INTENT(IN)           :: src_owner(:)
    INTEGER, INTENT(IN)           :: src_global_index(:)
    LOGICAL, OPTIONAL, INTENT(IN) :: inplace
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE interface_setup_comm_pattern

  SUBROUTINE interface_setup_comm_pattern2(p_pat, comm, recv_msg, send_msg, &
       glb2loc_index_recv, glb2loc_index_send, inplace)
    IMPORT t_comm_pattern, xfer_list, t_glb2loc_index_lookup
    CLASS(t_comm_pattern), TARGET, INTENT(OUT) :: p_pat
    INTEGER, INTENT(in) :: comm
    TYPE(xfer_list), INTENT(in), CONTIGUOUS :: recv_msg(:), send_msg(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index_recv, &
         glb2loc_index_send
    LOGICAL, OPTIONAL, INTENT(in) :: inplace
  END SUBROUTINE interface_setup_comm_pattern2

  SUBROUTINE interface_delete_comm_pattern(p_pat)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT) :: p_pat
  END SUBROUTINE interface_delete_comm_pattern

  SUBROUTINE interface_exchange_data_r3d(p_pat, lacc, recv, send, add)
    IMPORT t_comm_pattern, dp
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT)   :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
  END SUBROUTINE interface_exchange_data_r3d

  SUBROUTINE interface_exchange_data_s3d(p_pat, lacc, recv, send, add)
    IMPORT t_comm_pattern, sp
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT)   :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
  END SUBROUTINE interface_exchange_data_s3d

  SUBROUTINE interface_exchange_data_i3d(p_pat, lacc, recv, send, add)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT)  :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    INTEGER, INTENT(INOUT), TARGET        :: recv(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)
  END SUBROUTINE interface_exchange_data_i3d

  SUBROUTINE interface_exchange_data_mult( &
    p_pat, lacc, ndim2tot, recv, send, nshift)
    IMPORT t_comm_pattern, t_ptr_3d
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT) :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    TYPE(t_ptr_3d), PTR_INTENT(IN) :: recv(:)
    TYPE(t_ptr_3d), PTR_INTENT(IN), OPTIONAL :: send(:)
    INTEGER, INTENT(IN)           :: ndim2tot
    INTEGER, OPTIONAL, INTENT(IN) :: nshift
  END SUBROUTINE interface_exchange_data_mult

  SUBROUTINE interface_exchange_data_mult_mixprec( &
    p_pat, lacc, nfields_dp, ndim2tot_dp, nfields_sp, ndim2tot_sp, recv_dp, send_dp, &
    recv_sp, send_sp, nshift)
    IMPORT t_comm_pattern, t_ptr_3d, t_ptr_3d_sp
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT) :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    TYPE(t_ptr_3d), PTR_INTENT(in), OPTIONAL :: recv_dp(:)
    TYPE(t_ptr_3d), PTR_INTENT(in), OPTIONAL :: send_dp(:)
    TYPE(t_ptr_3d_sp), PTR_INTENT(in), OPTIONAL :: recv_sp(:)
    TYPE(t_ptr_3d_sp), PTR_INTENT(in), OPTIONAL :: send_sp(:)
    INTEGER, INTENT(IN)           :: nfields_dp, ndim2tot_dp, nfields_sp, &
         ndim2tot_sp
    INTEGER, OPTIONAL, INTENT(IN) :: nshift
  END SUBROUTINE interface_exchange_data_mult_mixprec

  SUBROUTINE interface_exchange_data_4de1(p_pat, lacc, nfields, ndim2tot, recv, send)
    IMPORT t_comm_pattern, dp
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT) :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT)           :: recv(:,:,:,:)
    REAL(dp), INTENT(IN   ), OPTIONAL :: send(:,:,:,:)
    INTEGER, INTENT(IN)           :: nfields, ndim2tot
  END SUBROUTINE interface_exchange_data_4de1

  SUBROUTINE interface_exchange_data_r2d(p_pat, lacc, recv, send, add, l_recv_exists)
    IMPORT t_comm_pattern, dp
    CLASS(t_comm_pattern), INTENT(INOUT), TARGET :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists
  END SUBROUTINE interface_exchange_data_r2d

  SUBROUTINE interface_exchange_data_s2d(p_pat, lacc, recv, send, add, l_recv_exists)
    IMPORT t_comm_pattern, sp
    CLASS(t_comm_pattern), INTENT(INOUT), TARGET :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists
  END SUBROUTINE interface_exchange_data_s2d

  SUBROUTINE interface_exchange_data_i2d(p_pat, lacc, recv, send, add, l_recv_exists)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), INTENT(INOUT), TARGET :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    INTEGER, INTENT(INOUT), TARGET        :: recv(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists
  END SUBROUTINE interface_exchange_data_i2d

  SUBROUTINE interface_exchange_data_l2d(p_pat, lacc, recv, send, l_recv_exists)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), INTENT(INOUT), TARGET :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    LOGICAL, OPTIONAL :: l_recv_exists
  END SUBROUTINE interface_exchange_data_l2d

  SUBROUTINE interface_exchange_data_l3d(p_pat, lacc, recv, send)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), TARGET, INTENT(INOUT)  :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
  END SUBROUTINE interface_exchange_data_l3d

  FUNCTION interface_get_np_recv(comm_pat)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), INTENT(IN) :: comm_pat
    INTEGER :: interface_get_np_recv
  END FUNCTION interface_get_np_recv

  FUNCTION interface_get_np_send(comm_pat)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), INTENT(IN) :: comm_pat
    INTEGER :: interface_get_np_send
  END FUNCTION interface_get_np_send

  SUBROUTINE interface_get_pelist_recv(comm_pat, pelist_recv)
    IMPORT t_comm_pattern
    CLASS(t_comm_pattern), INTENT(IN) :: comm_pat
    INTEGER, INTENT(OUT) :: pelist_recv(:)
  END SUBROUTINE interface_get_pelist_recv

  SUBROUTINE interface_delete_comm_pattern_collection(pattern_collection)
    IMPORT t_comm_pattern_collection
    CLASS(t_comm_pattern_collection), INTENT(INOUT) :: pattern_collection
  END SUBROUTINE interface_delete_comm_pattern_collection

  SUBROUTINE interface_exchange_data_grf(p_pat_coll, lacc, nfields, ndim2tot, &
       recv, send)
    IMPORT t_comm_pattern_collection, t_ptr_3d
    CLASS(t_comm_pattern_collection), TARGET, INTENT(INOUT) :: p_pat_coll
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    INTEGER, INTENT(IN) :: nfields
    INTEGER, INTENT(IN) :: ndim2tot
    ! recv itself is intent(in), but the pointed to data will be modified
    TYPE(t_ptr_3d), PTR_INTENT(in) :: recv(nfields), send(nfields)
  END SUBROUTINE interface_exchange_data_grf
END INTERFACE

END MODULE mo_communication_types
