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

! This module provides the communication routines
! for parallel runs

! (GZ, 2013-08-30): So far, the Cray compiler is the only one for which an OpenMP parallelization
! of copying data into / back from the MPI-buffer seems to give a benefit. Further compilers may
! be added here once the OpenMP implementation is sufficiently efficient
#if ((defined(_CRAYFTN) && !defined(_OPENACC)) || defined(__INTEL_COMPILER))
#define __OMPPAR_COPY__
#endif

!----------------------------
#include "icon_definitions.inc"
!----------------------------
MODULE mo_communication
!-------------------------------------------------------------------------
!
!    ProTeX FORTRAN source: Style 2
!    modified for ICON project, DWD/MPI-M 2006
!
!-------------------------------------------------------------------------
!
!
!
USE mo_impl_constants,       ONLY: SUCCESS
USE mo_scatter_pattern_base, ONLY: t_ScatterPattern, t_ScatterPatternPtr, deleteScatterPattern
USE mo_kind,                 ONLY: dp, sp
USE mo_exception,            ONLY: finish
USE mo_mpi,                  ONLY: p_send, p_recv, p_barrier, &
     &                             p_comm_work, p_pe_work, p_n_work, p_gatherv, &
     &                             p_alltoallv, p_alltoall, process_mpi_root_id, &
     &                             p_bcast, p_comm_is_intercomm, &
     &                             p_comm_remote_size, p_allgather, &
     &                             p_allgatherv, MPI_COMM_NULL
USE mo_parallel_config,      ONLY: blk_no, idx_no, idx_1d, nproma
USE mo_util_sort,            ONLY: quicksort
#ifdef __SX__
USE mo_util_sort,            ONLY: radixsort
#endif
USE mo_util_string,          ONLY: int2string
USE mo_fortran_tools,        ONLY: t_ptr_3d, t_ptr_3d_sp
USE mo_communication_types,  ONLY: t_comm_pattern, t_comm_pattern_collection, &
  &                                t_p_comm_pattern
USE mo_communication_types,  ONLY: t_comm_pattern, t_comm_pattern_collection, &
  &                                t_p_comm_pattern, xfer_list


IMPLICIT NONE

PRIVATE

!modules interface-------------------------------------------
!subroutines
PUBLIC :: blk_no, idx_no, idx_1d
PUBLIC :: delete_comm_pattern, exchange_data,                      &
          exchange_data_mult, exchange_data_grf,                   &
          exchange_data_4de1, exchange_data_mult_mixprec,          &
          get_np_recv, get_np_send, get_pelist_recv
PUBLIC :: t_comm_pattern
PUBLIC :: t_p_comm_pattern

PUBLIC :: t_comm_gather_pattern
PUBLIC :: setup_comm_gather_pattern
PUBLIC :: delete_comm_gather_pattern

PUBLIC :: t_comm_allgather_pattern
PUBLIC :: setup_comm_allgather_pattern
PUBLIC :: delete_comm_allgather_pattern

PUBLIC :: t_ScatterPattern, t_ScatterPatternPtr, makeScatterPattern, deleteScatterPattern

PUBLIC :: t_comm_pattern_collection, delete_comm_pattern_collection

PUBLIC :: ASSIGNMENT(=)

PUBLIC :: xfer_list

!
!------------------------------------------------------------------------------------------------
!
! MoHa: When changing t_comm_gather_pattern remember to adjust
!       copy_t_comm_gather_pattern accordingly!
TYPE t_comm_gather_pattern

  PRIVATE

  INTEGER, ALLOCATABLE :: collector_pes(:) ! ranks of collector processes
  INTEGER, ALLOCATABLE :: collector_size(:) ! total number of points per
                                            ! collector
  INTEGER, ALLOCATABLE :: collector_send_size(:) ! local number of points per
                                                 ! collector
  INTEGER, ALLOCATABLE :: loc_index(:) ! local indices of all points that
                                       ! need to be sent to a collector
  INTEGER, ALLOCATABLE :: recv_buffer_reorder(:) ! once the data is received
                                                 ! on the collectors, it has
                                                 ! to be reordered according
                                                 ! to this array
  INTEGER, ALLOCATABLE :: recv_buffer_reorder_fill(:) ! once the data is received
                                                      ! on the collectors, it has
                                                      ! to be reordered according
                                                      ! to this array (it leaves
                                                      ! holes for missing values)

  INTEGER, ALLOCATABLE :: recv_pes(:) ! ranks from which data is to be received
  INTEGER, ALLOCATABLE :: recv_size(:) ! number of remote points received
  INTEGER :: global_size ! global size of the array that is to be gathered
END TYPE t_comm_gather_pattern

INTERFACE ASSIGNMENT(=)
  MODULE PROCEDURE copy_t_comm_gather_pattern
END INTERFACE

!
!------------------------------------------------------------------------------------------------
!

TYPE t_comm_allgather_pattern
  TYPE(t_comm_gather_pattern), POINTER :: gather_pattern
  INTEGER :: intercomm
END TYPE t_comm_allgather_pattern


!--------------------------------------------------------------------------------------------------
!

INTERFACE exchange_data
   MODULE PROCEDURE exchange_data_r3d
   MODULE PROCEDURE exchange_data_s3d
   MODULE PROCEDURE exchange_data_i3d
   MODULE PROCEDURE exchange_data_r2d
   MODULE PROCEDURE exchange_data_s2d
   MODULE PROCEDURE exchange_data_i2d
   MODULE PROCEDURE exchange_data_l2d
   MODULE PROCEDURE exchange_data_l3d
   MODULE PROCEDURE gather_r_2d_deblock
   MODULE PROCEDURE gather_r_1d_deblock
   MODULE PROCEDURE gather_s_1d_deblock
   MODULE PROCEDURE gather_i_2d_deblock
   MODULE PROCEDURE gather_i_1d_deblock
   MODULE PROCEDURE allgather_r_1d_deblock
   MODULE PROCEDURE allgather_i_1d_deblock
END INTERFACE

INTERFACE two_phase_gather_first
  MODULE PROCEDURE two_phase_gather_first_r
  MODULE PROCEDURE two_phase_gather_first_s
  MODULE PROCEDURE two_phase_gather_first_i
END INTERFACE

INTERFACE two_phase_gather_second
  MODULE PROCEDURE two_phase_gather_second_r
  MODULE PROCEDURE two_phase_gather_second_s
  MODULE PROCEDURE two_phase_gather_second_i
END INTERFACE

CHARACTER(*), PARAMETER :: modname = "mo_communication"

!-------------------------------------------------------------------------

CONTAINS

  !-------------------------------------------------------------------------
  !
  !

  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------


  SUBROUTINE setup_comm_gather_pattern(global_size, owner_local, glb_index, &
    &                                  gather_pattern, disable_consistency_check)
    INTEGER, INTENT(IN) :: global_size, owner_local(:), glb_index(:)
    TYPE(t_comm_gather_pattern), INTENT(INOUT) :: gather_pattern
    LOGICAL, INTENT(IN), OPTIONAL :: disable_consistency_check

    INTEGER :: num_collectors
    LOGICAL, ALLOCATABLE :: pack_mask(:)
    INTEGER :: num_local_points, num_points_per_coll
    INTEGER, ALLOCATABLE :: packed_glb_index(:)
    INTEGER :: coll_stride
    INTEGER :: num_send_per_process(p_n_work), num_recv_per_process(p_n_work)
    INTEGER :: send_displ(p_n_work+1), recv_displ(p_n_work)
    INTEGER, ALLOCATABLE :: send_buffer(:), recv_buffer(:)
    INTEGER :: num_recv, num_recv_points
    INTEGER :: i, j, n
    LOGICAL :: consistency_check
    CHARACTER(*), PARAMETER :: routine = modname//":setup_comm_gather_pattern"

    gather_pattern%global_size = global_size

    ! determine collector ranks and the data associated to each collector
    num_collectors = NINT(SQRT(REAL(p_n_work)))
    IF (ALLOCATED(gather_pattern%collector_pes)) &
      DEALLOCATE(gather_pattern%collector_pes)
    IF (ALLOCATED(gather_pattern%collector_size)) &
      DEALLOCATE(gather_pattern%collector_size)
    IF (ALLOCATED(gather_pattern%collector_send_size)) &
      DEALLOCATE(gather_pattern%collector_send_size)
    ALLOCATE(gather_pattern%collector_pes(num_collectors), &
      &      gather_pattern%collector_size(num_collectors), &
      &      gather_pattern%collector_send_size(num_collectors))
    coll_stride = (p_n_work + num_collectors - 1) / &
      &           num_collectors
    num_points_per_coll = (global_size + num_collectors - 1) / num_collectors
    DO i = 1, num_collectors
      ! set collector ranks
      gather_pattern%collector_pes(i) = (i-1) * coll_stride
    END DO

    ! mask for all locally owned points
    ALLOCATE(pack_mask(SIZE(owner_local)))
    pack_mask(:) = owner_local(:) == p_pe_work
    num_local_points = COUNT(pack_mask(:))

    ! determine local indices of all points that need to be sent to a collector
    IF (ALLOCATED(gather_pattern%loc_index)) &
      DEALLOCATE(gather_pattern%loc_index)
    ALLOCATE(gather_pattern%loc_index(num_local_points), &
      &      packed_glb_index(num_local_points))
    packed_glb_index(:) = PACK(glb_index(:), pack_mask(:))
    gather_pattern%loc_index(:) = PACK((/(i, i = 1, &
      &                                SIZE(owner_local))/), pack_mask(:))

    DEALLOCATE(pack_mask)

    ! sort loc_index according to the respective global indices
#ifdef __SX__
    CALL radixsort(packed_glb_index(:), gather_pattern%loc_index(:))
#else
    CALL quicksort(packed_glb_index(:), gather_pattern%loc_index(:))
#endif
    ! determine number of points that need to be sent to each collector
    gather_pattern%collector_send_size(:) = 0
    DO i = 1, num_local_points
      n = 1 + (packed_glb_index(i) - 1) / num_points_per_coll
      gather_pattern%collector_send_size(n) = &
        gather_pattern%collector_send_size(n) + 1
    END DO

    ! generate send and receive counts for all processes
    num_send_per_process(:) = 0
    DO i = 1, num_collectors
      num_send_per_process(gather_pattern%collector_pes(i)+1) = &
        gather_pattern%collector_send_size(i)
    END DO
    CALL p_alltoall(num_send_per_process(:), num_recv_per_process(:), &
      &             p_comm_work)
    num_recv_points = SUM(num_recv_per_process(:))

    ! exchange number of points per collector
    IF (p_pe_work == gather_pattern%collector_pes(1)) THEN
      gather_pattern%collector_size(1) = num_recv_points
      DO i = 2, num_collectors
        CALL p_recv(gather_pattern%collector_size(i), &
          &         gather_pattern%collector_pes(i), 0, 1, p_comm_work)
      END DO
    ELSE IF (ANY(gather_pattern%collector_pes(:) == p_pe_work)) THEN
      CALL p_send(num_recv_points, gather_pattern%collector_pes(1), 0, 1, &
        &         p_comm_work)
    END IF
    CALL p_bcast(gather_pattern%collector_size(:), &
      &          gather_pattern%collector_pes(1), p_comm_work)

    ! number of messages to be received (is 0 on non-collector processes)
    num_recv = COUNT(num_recv_per_process(:) /= 0)
    IF (ALLOCATED(gather_pattern%recv_pes)) &
      DEALLOCATE(gather_pattern%recv_pes)
    IF (ALLOCATED(gather_pattern%recv_size)) &
      DEALLOCATE(gather_pattern%recv_size)
    ALLOCATE(gather_pattern%recv_pes(num_recv), &
      &      gather_pattern%recv_size(num_recv))
    num_recv = 0
    DO i = 1, p_n_work
      IF (num_recv_per_process(i) /= 0) THEN
        num_recv = num_recv + 1
        gather_pattern%recv_pes(num_recv) = i - 1
        gather_pattern%recv_size(num_recv) = num_recv_per_process(i)
      END IF
    END DO

    ! generate send and receive displacement and fill send buffer
    ! remark: at first the content of send_displ is shifted by one element
    !         to the back, this is done in order to ease the copying of data
    !         into the send buffer
    send_displ(1:2) = 0
    recv_displ(1) = 0
    DO i = 2, p_n_work
      send_displ(i+1) = send_displ(i) + num_send_per_process(i-1)
      recv_displ(i)   = recv_displ(i-1) + num_recv_per_process(i-1)
    END DO
    ALLOCATE(send_buffer(SUM(num_send_per_process(:))), &
      &      recv_buffer(num_recv_points))

    DO i = 1, num_local_points
      j = 2 + gather_pattern%collector_pes(1 + (packed_glb_index(i)-1) / &
        &                                  num_points_per_coll)
      send_displ(j) = send_displ(j) + 1
      send_buffer(send_displ(j)) = packed_glb_index(i)
    END DO

    DEALLOCATE(packed_glb_index)

    ! collect the global indices from all processes
    CALL p_alltoallv(send_buffer, num_send_per_process, send_displ, &
      &              recv_buffer, num_recv_per_process, recv_displ, &
      &              p_comm_work)

    ! compute the final position of received data on the collectors
    IF (ALLOCATED(gather_pattern%recv_buffer_reorder)) &
      DEALLOCATE(gather_pattern%recv_buffer_reorder)
    ALLOCATE(gather_pattern%recv_buffer_reorder(num_recv_points))
    gather_pattern%recv_buffer_reorder(:) = (/(i, i = 1, num_recv_points)/)
#ifdef __SX__
    CALL radixsort(recv_buffer(:), gather_pattern%recv_buffer_reorder(:))
#else
    CALL quicksort(recv_buffer(:), gather_pattern%recv_buffer_reorder(:))
#endif
    IF (ALLOCATED(gather_pattern%recv_buffer_reorder_fill)) &
      DEALLOCATE(gather_pattern%recv_buffer_reorder_fill)
    ALLOCATE(gather_pattern%recv_buffer_reorder_fill(num_recv_points))
    gather_pattern%recv_buffer_reorder_fill(:) = &
      MOD(recv_buffer-1, num_points_per_coll)+1

    IF (PRESENT(disable_consistency_check)) THEN
      consistency_check = .NOT. disable_consistency_check
    ELSE
      consistency_check = .TRUE.
    END IF

    ! consistency check
    IF (consistency_check) THEN

      ! check whether a point is owned by multiple processes
      DO i = 2, num_recv_points
        IF (recv_buffer(i) == recv_buffer(i-1)) &
          CALL finish(routine, &
          &         "One or more points are owned by multiple processes")
      END DO

      ! check if all points have an owner
      IF (SUM(gather_pattern%collector_size(:)) /= global_size) &
        CALL finish(routine, int2string(SUM(gather_pattern%collector_size(:)))//" out of "//int2string(global_size)//&
        &" points have no owner!")
    END IF

    DEALLOCATE(send_buffer, recv_buffer)

  END SUBROUTINE setup_comm_gather_pattern


  !-------------------------------------------------------------------------


  SUBROUTINE setup_comm_allgather_pattern(gather_pattern, intercomm, &
    &                                     allgather_pattern)
    TYPE(t_comm_gather_pattern), TARGET, INTENT(INOUT) :: gather_pattern
    INTEGER, OPTIONAL, INTENT(IN) :: intercomm
    TYPE(t_comm_allgather_pattern), INTENT(INOUT) :: allgather_pattern

    IF (PRESENT(intercomm)) THEN
      ! check whether intercomm is really a intercommunicator
      IF (.NOT. p_comm_is_intercomm(intercomm)) &
        CALL finish("setup_comm_allgather_pattern", "invalid intercomm")
      allgather_pattern%intercomm = intercomm
    ELSE
      allgather_pattern%intercomm = MPI_COMM_NULL
    END IF
    allgather_pattern%gather_pattern => gather_pattern
  END SUBROUTINE setup_comm_allgather_pattern


  !-------------------------------------------------------------------------
  !
  !! Deletes a communication pattern, i.e. deallocates all arrays
  !! and sets all other members to 0
  !
  SUBROUTINE delete_comm_pattern(p_pat)

    CLASS(t_comm_pattern), POINTER :: p_pat

    CALL p_pat%delete()

    DEALLOCATE(p_pat)

  END SUBROUTINE delete_comm_pattern


  !-------------------------------------------------------------------------


  SUBROUTINE delete_comm_pattern_collection(pattern_collection)

    CLASS(t_comm_pattern_collection), POINTER :: pattern_collection

    CALL pattern_collection%delete()

    DEALLOCATE(pattern_collection)

  END SUBROUTINE delete_comm_pattern_collection


  !-------------------------------------------------------------------------


  SUBROUTINE delete_comm_gather_pattern(gather_pattern)
    TYPE(t_comm_gather_pattern), INTENT(INOUT) :: gather_pattern

    IF (ALLOCATED(gather_pattern%collector_pes)) &
      DEALLOCATE(gather_pattern%collector_pes)
    IF (ALLOCATED(gather_pattern%collector_size)) &
      DEALLOCATE(gather_pattern%collector_size)
    IF (ALLOCATED(gather_pattern%collector_send_size)) &
      DEALLOCATE(gather_pattern%collector_send_size)
    IF (ALLOCATED(gather_pattern%loc_index)) &
      DEALLOCATE(gather_pattern%loc_index)
    IF (ALLOCATED(gather_pattern%recv_buffer_reorder)) &
      DEALLOCATE(gather_pattern%recv_buffer_reorder)
    IF (ALLOCATED(gather_pattern%recv_buffer_reorder_fill)) &
      DEALLOCATE(gather_pattern%recv_buffer_reorder_fill)
    IF (ALLOCATED(gather_pattern%recv_pes)) &
      DEALLOCATE(gather_pattern%recv_pes)
    IF (ALLOCATED(gather_pattern%recv_size)) &
      DEALLOCATE(gather_pattern%recv_size)
  END SUBROUTINE delete_comm_gather_pattern


  !-------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE copy_t_comm_gather_pattern(out_arg, in_arg)

    TYPE(t_comm_gather_pattern), INTENT(OUT) :: out_arg
    TYPE(t_comm_gather_pattern), INTENT(IN) :: in_arg

    IF (ALLOCATED(in_arg%collector_pes)) THEN
      ALLOCATE(out_arg%collector_pes(SIZE(in_arg%collector_pes)))
      out_arg%collector_pes(:) = in_arg%collector_pes(:)
    END IF
    IF (ALLOCATED(in_arg%collector_size)) THEN
      ALLOCATE(out_arg%collector_size(SIZE(in_arg%collector_size)))
      out_arg%collector_size(:) = in_arg%collector_size(:)
    END IF
    IF (ALLOCATED(in_arg%collector_send_size)) THEN
      ALLOCATE(out_arg%collector_send_size(SIZE(in_arg%collector_send_size)))
      out_arg%collector_send_size(:) = in_arg%collector_send_size(:)
    END IF
    IF (ALLOCATED(in_arg%loc_index)) THEN
      ALLOCATE(out_arg%loc_index(SIZE(in_arg%loc_index)))
      out_arg%loc_index(:) = in_arg%loc_index(:)
    END IF
    IF (ALLOCATED(in_arg%recv_buffer_reorder)) THEN
      ALLOCATE(out_arg%recv_buffer_reorder(SIZE(in_arg%recv_buffer_reorder)))
      out_arg%recv_buffer_reorder(:) = in_arg%recv_buffer_reorder(:)
    END IF
    IF (ALLOCATED(in_arg%recv_buffer_reorder_fill)) THEN
      ALLOCATE( &
           out_arg%recv_buffer_reorder_fill(SIZE(in_arg%recv_buffer_reorder_fill)))
      out_arg%recv_buffer_reorder_fill(:) = in_arg%recv_buffer_reorder_fill(:)
    END IF
    IF (ALLOCATED(in_arg%recv_pes)) THEN
      ALLOCATE(out_arg%recv_pes(SIZE(in_arg%recv_pes)))
      out_arg%recv_pes(:) = in_arg%recv_pes(:)
    END IF
    IF (ALLOCATED(in_arg%recv_size)) THEN
      ALLOCATE(out_arg%recv_size(SIZE(in_arg%recv_size)))
      out_arg%recv_size(:) = in_arg%recv_size(:)
    END IF
    out_arg%global_size = in_arg%global_size

  END SUBROUTINE copy_t_comm_gather_pattern


  !-------------------------------------------------------------------------

  SUBROUTINE delete_comm_allgather_pattern(allgather_pattern)
    TYPE(t_comm_allgather_pattern), INTENT(INOUT) :: allgather_pattern

    RETURN

  END SUBROUTINE delete_comm_allgather_pattern


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_r3d(p_pat, lacc, recv, send, add)

    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)

    CALL p_pat%exchange_data_r3d(lacc, recv, send, add)

  END SUBROUTINE exchange_data_r3d


  SUBROUTINE exchange_data_s3d(p_pat, lacc, recv, send, add)

    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)

    CALL p_pat%exchange_data_s3d(lacc, recv, send, add)

  END SUBROUTINE exchange_data_s3d


  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_i3d(p_pat, lacc, recv, send, add)

    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    INTEGER, INTENT(INOUT), TARGET        :: recv(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:,:)

    CALL p_pat%exchange_data_i3d(lacc, recv, send, add)

  END SUBROUTINE exchange_data_i3d

  SUBROUTINE exchange_data_l3d(p_pat, lacc, recv, send)

    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:,:)

    CALL p_pat%exchange_data_l3d(lacc, recv, send)

  END SUBROUTINE exchange_data_l3d

  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  SUBROUTINE exchange_data_4de1(p_pat, lacc, nfields, ndim2tot, recv, send)

    CLASS(t_comm_pattern), POINTER :: p_pat

    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT)           :: recv(:,:,:,:)
    REAL(dp), INTENT(IN   ), OPTIONAL :: send(:,:,:,:)
    INTEGER, INTENT(IN)           :: nfields, ndim2tot

    CALL p_pat%exchange_data_4de1(lacc, nfields, ndim2tot, recv, send)

  END SUBROUTINE exchange_data_4de1


  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  SUBROUTINE exchange_data_grf(p_pat_coll, lacc, nfields, ndim2tot, recv1, send1, &
    recv2, send2, recv3, send3, recv4, send4, &
    recv5, send5, recv6, send6, recv4d1, send4d1, &
    recv4d2, send4d2)

    CLASS(t_comm_pattern_collection), POINTER :: p_pat_coll

    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4d1(:,:,:,:), &
      recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), recv4d2(:,:,:,:)
    ! Note: the last index of the send fields corresponds to the dimension of p_pat
    ! On the other hand, they are not blocked and have the vertical index first
    REAL(dp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1(:,:,:), send2(:,:,:), send3(:,:,:), send4d1(:,:,:,:), &
      send4(:,:,:), send5(:,:,:), send6(:,:,:), send4d2(:,:,:,:)

    INTEGER, INTENT(IN)           :: nfields  ! total number of input fields
    INTEGER, INTENT(IN)           :: ndim2tot ! sum of vertical levels of input fields

    TYPE(t_ptr_3d) :: recv(nfields), send(nfields)

    INTEGER :: n, n4d

    ! Set pointers to input fields
    IF (PRESENT(recv4d1) .AND. .NOT. PRESENT(recv4d2)) THEN
      DO n = 1, nfields
        recv(n)%p => recv4d1(:,:,:,n)
        send(n)%p => send4d1(:,:,:,n)
      ENDDO
    ELSE IF (PRESENT(recv4d1) .AND. PRESENT(recv4d2)) THEN
      n4d = nfields/2
      DO n = 1, n4d
        recv(n)%p => recv4d1(:,:,:,n)
        send(n)%p => send4d1(:,:,:,n)
      ENDDO
      DO n = 1, n4d
        recv(n4d+n)%p => recv4d2(:,:,:,n)
        send(n4d+n)%p => send4d2(:,:,:,n)
      ENDDO
    ELSE
      IF (PRESENT(recv1)) THEN
        recv(1)%p => recv1
        send(1)%p => send1
      ENDIF
      IF (PRESENT(recv2)) THEN
        recv(2)%p => recv2
        send(2)%p => send2
      ENDIF
      IF (PRESENT(recv3)) THEN
        recv(3)%p => recv3
        send(3)%p => send3
      ENDIF
      IF (PRESENT(recv4)) THEN
        recv(4)%p => recv4
        send(4)%p => send4
      ENDIF
      IF (PRESENT(recv5)) THEN
        recv(5)%p => recv5
        send(5)%p => send5
      ENDIF
      IF (PRESENT(recv6)) THEN
        recv(6)%p => recv6
        send(6)%p => send6
      ENDIF
    ENDIF

    CALL p_pat_coll%exchange_data_grf(lacc, nfields, ndim2tot, recv, send)

  END SUBROUTINE exchange_data_grf


  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Interface for 2D arrays for exchange_data.
  !!
  !! Just reshapes the arrays and calls exchange_data.
  !!
  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_r2d(p_pat, lacc, recv, send, add, l_recv_exists)
    !
    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(dp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CALL p_pat%exchange_data_r2d(lacc, recv, send, add, l_recv_exists)

  END SUBROUTINE exchange_data_r2d

  SUBROUTINE exchange_data_s2d(p_pat, lacc, recv, send, add, l_recv_exists)
    !
    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(sp), INTENT(INOUT), TARGET        :: recv(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    REAL(sp), INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CALL p_pat%exchange_data_s2d(lacc, recv, send, add, l_recv_exists)

  END SUBROUTINE exchange_data_s2d


  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE exchange_data_i2d(p_pat, lacc, recv, send, add, l_recv_exists)
    !
    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    INTEGER, INTENT(INOUT), TARGET        :: recv(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    INTEGER, INTENT(IN), OPTIONAL, TARGET :: add (:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CALL p_pat%exchange_data_i2d(lacc, recv, send, add, l_recv_exists)

  END SUBROUTINE exchange_data_i2d

  SUBROUTINE exchange_data_l2d(p_pat, lacc, recv, send, l_recv_exists)
    !
    CLASS(t_comm_pattern), POINTER :: p_pat
    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    LOGICAL, INTENT(INOUT), TARGET        :: recv(:,:)
    LOGICAL, INTENT(IN), OPTIONAL, TARGET :: send(:,:)
    LOGICAL, OPTIONAL :: l_recv_exists

    CALL p_pat%exchange_data_l2d(lacc, recv, send, l_recv_exists)

  END SUBROUTINE exchange_data_l2d

  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  SUBROUTINE exchange_data_mult(p_pat, lacc, nfields, ndim2tot, recv1, send1, recv2, send2,   &
    recv3, send3, recv4, send4,  recv5, send5, recv6, send6, recv7, send7,              &
    recv4d, send4d, nshift, recv3d_arr, send3d_arr)

    CLASS(t_comm_pattern), POINTER, INTENT(in) :: p_pat

    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1(:,:,:), recv2(:,:,:), recv3(:,:,:), recv4(:,:,:), recv5(:,:,:), recv6(:,:,:), &
      recv7(:,:,:), recv4d(:,:,:,:)
    REAL(dp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1(:,:,:), send2(:,:,:), send3(:,:,:), send4(:,:,:), send5(:,:,:), send6(:,:,:), &
      send7(:,:,:), send4d(:,:,:,:)

    INTEGER, INTENT(IN)           :: nfields, ndim2tot
    TYPE(t_ptr_3d), INTENT(   IN), TARGET, OPTIONAL :: recv3d_arr(:)
    TYPE(t_ptr_3d), INTENT(INOUT), TARGET, OPTIONAL :: send3d_arr(:)
    INTEGER, OPTIONAL, INTENT(IN) :: nshift

    CHARACTER(len=*), PARAMETER :: routine = modname//"::exchange_data_mult"
    TYPE(t_ptr_3d) :: recv(nfields), send(nfields)

    INTEGER :: i, nf4d
    LOGICAL :: lsend

    !-----------------------------------------------------------------------
    lsend     = .FALSE.

    ! Set pointers to input fields
    IF (PRESENT(recv4d)) THEN
      nf4d = SIZE(recv4d,4)
      DO i = 1, nf4d
        recv(i)%p => recv4d(:,:,:,i)
      ENDDO
      IF (PRESENT(send4d)) THEN ! all 4D fields must have the same dimensions
        DO i = 1, nf4d
          send(i)%p => send4d(:,:,:,i)
        ENDDO
        lsend = .TRUE.
      ENDIF
    ELSE
      nf4d = 0
    ENDIF


    ! Set pointers to input fields
    IF (PRESENT(recv3d_arr)) THEN
      DO i = 1, SIZE(recv3d_arr)
        recv(i+nf4d)%p => recv3d_arr(i)%p
      ENDDO
      IF (PRESENT(send3d_arr)) THEN
        DO i = 1, SIZE(recv3d_arr)
          send(i+nf4d)%p => send3d_arr(i)%p
        ENDDO
        lsend = .TRUE.
      ENDIF
      nf4d = nf4d + SIZE(recv3d_arr)
    ENDIF


    IF (PRESENT(recv1)) THEN
      recv(nf4d+1)%p => recv1
      IF (PRESENT(send1)) THEN
        send(nf4d+1)%p => send1
        lsend = .TRUE.
      ENDIF
      IF (PRESENT(recv2)) THEN
        recv(nf4d+2)%p => recv2
        IF (lsend) send(nf4d+2)%p => send2
        IF (PRESENT(recv3)) THEN
          recv(nf4d+3)%p => recv3
          IF (lsend) send(nf4d+3)%p => send3
          IF (PRESENT(recv4)) THEN
            recv(nf4d+4)%p => recv4
            IF (lsend) send(nf4d+4)%p => send4
            IF (PRESENT(recv5)) THEN
              recv(nf4d+5)%p => recv5
              IF (lsend) send(nf4d+5)%p => send5
              IF (PRESENT(recv6)) THEN
                recv(nf4d+6)%p => recv6
                IF (lsend) send(nf4d+6)%p => send6
                IF (PRESENT(recv7)) THEN
                  recv(nf4d+7)%p => recv7
                  IF (lsend) send(nf4d+7)%p => send7
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF

    IF (lsend) THEN
      CALL p_pat%exchange_data_mult(lacc, ndim2tot, recv, send, nshift)
    ELSE
      CALL p_pat%exchange_data_mult(lacc, ndim2tot, recv, nshift=nshift)
    END IF

  END SUBROUTINE exchange_data_mult



  !! Does data exchange according to a communication pattern (in p_pat).
  !!
  SUBROUTINE exchange_data_mult_mixprec(p_pat, lacc, nfields_dp, ndim2tot_dp, nfields_sp, ndim2tot_sp,        &
    recv1_dp, send1_dp, recv2_dp, send2_dp, recv3_dp, send3_dp, recv4_dp, send4_dp, recv5_dp, send5_dp, &
    recv1_sp, send1_sp, recv2_sp, send2_sp, recv3_sp, send3_sp, recv4_sp, send4_sp, recv5_sp, send5_sp, &
    recv4d_dp, send4d_dp, recv4d_sp, send4d_sp, nshift)

    CLASS(t_comm_pattern), INTENT(INOUT) :: p_pat

    LOGICAL, INTENT(IN) :: lacc ! If true, use openacc
    REAL(dp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1_dp(:,:,:), recv2_dp(:,:,:), recv3_dp(:,:,:), recv4_dp(:,:,:), recv5_dp(:,:,:), recv4d_dp(:,:,:,:)
    REAL(dp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1_dp(:,:,:), send2_dp(:,:,:), send3_dp(:,:,:), send4_dp(:,:,:), send5_dp(:,:,:), send4d_dp(:,:,:,:)

    REAL(sp), INTENT(INOUT), TARGET, OPTIONAL ::  &
      recv1_sp(:,:,:), recv2_sp(:,:,:), recv3_sp(:,:,:), recv4_sp(:,:,:), recv5_sp(:,:,:), recv4d_sp(:,:,:,:)
    REAL(sp), INTENT(IN   ), TARGET, OPTIONAL ::  &
      send1_sp(:,:,:), send2_sp(:,:,:), send3_sp(:,:,:), send4_sp(:,:,:), send5_sp(:,:,:), send4d_sp(:,:,:,:)

    INTEGER, INTENT(IN)           :: nfields_dp, ndim2tot_dp, nfields_sp, ndim2tot_sp
    INTEGER, OPTIONAL, INTENT(IN) :: nshift


    TYPE(t_ptr_3d) :: recv_dp(nfields_dp), send_dp(nfields_dp)
    TYPE(t_ptr_3d_sp) :: recv_sp(nfields_sp), send_sp(nfields_sp)

    INTEGER :: i, nf4d_dp, nf4d_sp
    LOGICAL :: lsend

    !-----------------------------------------------------------------------

    lsend     = .FALSE.

    ! Set pointers to input fields
    IF (PRESENT(recv4d_dp)) THEN
      nf4d_dp = SIZE(recv4d_dp,4)
      DO i = 1, nf4d_dp
        recv_dp(i)%p => recv4d_dp(:,:,:,i)
      ENDDO
      IF (PRESENT(send4d_dp)) THEN ! all 4D fields must have the same dimensions
        DO i = 1, nf4d_dp
          send_dp(i)%p => send4d_dp(:,:,:,i)
        ENDDO
        lsend = .TRUE.
      ENDIF
    ELSE
      nf4d_dp = 0
    ENDIF
    IF (PRESENT(recv4d_sp)) THEN
      nf4d_sp = SIZE(recv4d_sp,4)
      DO i = 1, nf4d_sp
        recv_sp(i)%p => recv4d_sp(:,:,:,i)
      ENDDO
      IF (PRESENT(send4d_sp)) THEN ! all 4D fields must have the same dimensions
        DO i = 1, nf4d_sp
          send_sp(i)%p => send4d_sp(:,:,:,i)
        ENDDO
        lsend = .TRUE.
      ENDIF
    ELSE
      nf4d_sp = 0
    ENDIF

    IF (PRESENT(recv1_dp)) THEN
      recv_dp(nf4d_dp+1)%p => recv1_dp
      IF (PRESENT(send1_dp)) THEN
        send_dp(nf4d_dp+1)%p => send1_dp
        lsend = .TRUE.
      ENDIF
      IF (PRESENT(recv2_dp)) THEN
        recv_dp(nf4d_dp+2)%p => recv2_dp
        IF (lsend) send_dp(nf4d_dp+2)%p => send2_dp
        IF (PRESENT(recv3_dp)) THEN
          recv_dp(nf4d_dp+3)%p => recv3_dp
          IF (lsend) send_dp(nf4d_dp+3)%p => send3_dp
          IF (PRESENT(recv4_dp)) THEN
            recv_dp(nf4d_dp+4)%p => recv4_dp
            IF (lsend) send_dp(nf4d_dp+4)%p => send4_dp
            IF (PRESENT(recv5_dp)) THEN
              recv_dp(nf4d_dp+5)%p => recv5_dp
              IF (lsend) send_dp(nf4d_dp+5)%p => send5_dp
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    IF (PRESENT(recv1_sp)) THEN
      recv_sp(nf4d_sp+1)%p => recv1_sp
      IF (PRESENT(send1_sp)) THEN
        send_sp(nf4d_sp+1)%p => send1_sp
        lsend = .TRUE.
      ENDIF
      IF (PRESENT(recv2_sp)) THEN
        recv_sp(nf4d_sp+2)%p => recv2_sp
        IF (lsend) send_sp(nf4d_sp+2)%p => send2_sp
        IF (PRESENT(recv3_sp)) THEN
          recv_sp(nf4d_sp+3)%p => recv3_sp
          IF (lsend) send_sp(nf4d_sp+3)%p => send3_sp
          IF (PRESENT(recv4_sp)) THEN
            recv_sp(nf4d_sp+4)%p => recv4_sp
            IF (lsend) send_sp(nf4d_sp+4)%p => send4_sp
            IF (PRESENT(recv5_sp)) THEN
              recv_sp(nf4d_sp+5)%p => recv5_sp
              IF (lsend) send_sp(nf4d_sp+5)%p => send5_sp
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF

    IF (lsend) THEN
      CALL p_pat%exchange_data_mult_mixprec(lacc, nfields_dp, ndim2tot_dp, &
           nfields_sp, ndim2tot_sp, recv_dp, send_dp, recv_sp, send_sp, nshift)
    ELSE
      CALL p_pat%exchange_data_mult_mixprec(lacc, nfields_dp, ndim2tot_dp, &
           nfields_sp, ndim2tot_sp, recv_dp=recv_dp, recv_sp=recv_sp, &
           nshift=nshift)
    END IF



  END SUBROUTINE exchange_data_mult_mixprec


  !-------------------------------------------------------------------------


  FUNCTION get_np_recv(comm_pat)
    CLASS(t_comm_pattern), POINTER :: comm_pat
    INTEGER :: get_np_recv

    get_np_recv = comm_pat%get_np_recv()
  END FUNCTION get_np_recv


  !-------------------------------------------------------------------------


  FUNCTION get_np_send(comm_pat)
    CLASS(t_comm_pattern), POINTER :: comm_pat
    INTEGER :: get_np_send

    get_np_send =  comm_pat%get_np_send()
  END FUNCTION get_np_send


  !-------------------------------------------------------------------------


  SUBROUTINE get_pelist_recv(comm_pat, pelist_recv)
    CLASS(t_comm_pattern), POINTER :: comm_pat
    INTEGER, INTENT(OUT) :: pelist_recv(:)

    CALL comm_pat%get_pelist_recv(pelist_recv)
  END SUBROUTINE get_pelist_recv


  !-------------------------------------------------------------------------


  SUBROUTINE gather_r_1d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nblk)
    REAL(dp), INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    REAL(dp), INTENT(INOUT) :: out_array(:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    REAL(dp), ALLOCATABLE :: send_buffer(:,:)
    REAL(dp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)
    CALL out_array_to_2d(out_array, SIZE(out_array))

    DEALLOCATE(send_buffer)

  CONTAINS

    SUBROUTINE out_array_to_2d(out_array_2d, n)
      INTEGER, INTENT(IN) :: n
      REAL(dp), INTENT(INOUT) :: out_array_2d(1,n)

      CALL two_phase_gather_second(recv_buffer_r=out_array_2d, &
        fill_value=fill_value, &
        gather_pattern=gather_pattern, &
        collector_buffer_r=collector_buffer)
    END SUBROUTINE out_array_to_2d

  END SUBROUTINE gather_r_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_s_1d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nblk)
    REAL(sp), INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    REAL(sp), INTENT(INOUT) :: out_array(:)
    REAL(sp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    REAL(sp), ALLOCATABLE :: send_buffer(:,:)
    REAL(sp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)
    CALL out_array_to_2d(out_array, SIZE(out_array))

    DEALLOCATE(send_buffer)

  CONTAINS

    SUBROUTINE out_array_to_2d(out_array_2d, n)
      INTEGER, INTENT(IN) :: n
      REAL(sp), INTENT(INOUT) :: out_array_2d(1,n)

      CALL two_phase_gather_second(recv_buffer_r=out_array_2d, &
        fill_value=fill_value, &
        gather_pattern=gather_pattern, &
        collector_buffer_r=collector_buffer)
    END SUBROUTINE out_array_to_2d

  END SUBROUTINE gather_s_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_i_1d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nblk)
    INTEGER, INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    INTEGER, INTENT(INOUT) :: out_array(:)
    INTEGER, INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    INTEGER, ALLOCATABLE :: send_buffer(:,:)
    INTEGER, POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_i=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_i=collector_buffer)
    CALL out_array_to_2d(out_array, SIZE(out_array))

    DEALLOCATE(send_buffer)

  CONTAINS

    SUBROUTINE out_array_to_2d(out_array_2d, n)
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(INOUT) :: out_array_2d(1,n)

      CALL two_phase_gather_second(recv_buffer_i=out_array_2d, &
        fill_value=fill_value, &
        gather_pattern=gather_pattern, &
        collector_buffer_i=collector_buffer)
    END SUBROUTINE out_array_to_2d

  END SUBROUTINE gather_i_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_r_2d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nlev, nblk)
    REAL(dp), INTENT(IN) :: in_array(:,:,:)
    ! dimension (global length, nlev); only required on root
    REAL(dp), INTENT(INOUT) :: out_array(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    REAL(dp), ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
    REAL(dp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, nlev, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    nlev = SIZE(in_array, 2)

    IF (SIZE(in_array, 1) /= nproma) &
      CALL finish("gather_r_2d_deblock", &
      &         "size of first dimension of in_array is not nproma")

    IF (nlev /= SIZE(out_array, 2) .AND. p_pe_work == process_mpi_root_id) &
      CALL finish("gather_r_2d_deblock", &
      &         "second size of in_array and out_array are not the same")

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    IF (SIZE(in_array, 1) * SIZE(in_array, 3) < num_send_points) &
      CALL finish("gather_r_2d_deblock", "in_array is too small")

    ALLOCATE(send_buffer(nlev, num_send_points))
    IF (p_pe_work == process_mpi_root_id) THEN
      ALLOCATE(recv_buffer(nlev, MERGE(gather_pattern%global_size, &
        &                              SUM(gather_pattern%collector_size(:)), &
        &                              PRESENT(fill_value))))
    ELSE
      ALLOCATE(recv_buffer(0,0))
    END IF

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(:,i) = in_array(idx, :, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)
    CALL two_phase_gather_second(recv_buffer_r=recv_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_r=collector_buffer)

    IF (p_pe_work == process_mpi_root_id) &
      out_array(1:SIZE(recv_buffer, 2),1:SIZE(recv_buffer, 1)) = &
      TRANSPOSE(recv_buffer(:,:))

    DEALLOCATE(send_buffer,recv_buffer)

  END SUBROUTINE gather_r_2d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE gather_i_2d_deblock(in_array, out_array, fill_value, gather_pattern)
    ! dimension (nproma, nlev, nblk)
    INTEGER, INTENT(IN) :: in_array(:,:,:)
    ! dimension (global length, nlev); only required on root
    INTEGER, INTENT(INOUT) :: out_array(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern

    INTEGER, ALLOCATABLE :: send_buffer(:,:), recv_buffer(:,:)
    INTEGER, POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, nlev, idx, blk

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    nlev = SIZE(in_array, 2)

    IF (SIZE(in_array, 1) /= nproma) &
      CALL finish("gather_i_2d_deblock", &
      &         "size of first dimension of in_array is not nproma")

    IF (nlev /= SIZE(out_array, 2) .AND. p_pe_work == process_mpi_root_id) &
      CALL finish("gather_i_2d_deblock", &
      &         "second size of in_array and out_array are not the same")

    num_send_points = SUM(gather_pattern%collector_send_size(:))
    IF (SIZE(in_array, 1) * SIZE(in_array, 3) < num_send_points) &
      CALL finish("gather_i_2d_deblock", "in_array is too small")

    ALLOCATE(send_buffer(nlev, num_send_points))
    IF (p_pe_work == process_mpi_root_id) THEN
      ALLOCATE(recv_buffer(nlev, MERGE(gather_pattern%global_size, &
        &                              SUM(gather_pattern%collector_size(:)), &
        &                              PRESENT(fill_value))))
    ELSE
      ALLOCATE(recv_buffer(0,0))
    END IF

    DO i = 1, SIZE(gather_pattern%loc_index(:))
      idx = idx_no(gather_pattern%loc_index(i))
      blk = blk_no(gather_pattern%loc_index(i))
      send_buffer(:,i) = in_array(idx, :, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_i=send_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_i=collector_buffer)
    CALL two_phase_gather_second(recv_buffer_i=recv_buffer, fill_value=fill_value,&
      gather_pattern=gather_pattern, &
      collector_buffer_i=collector_buffer)

    IF (p_pe_work == process_mpi_root_id) &
      out_array(1:SIZE(recv_buffer, 2),1:SIZE(recv_buffer, 1)) = &
      TRANSPOSE(recv_buffer(:,:))

    DEALLOCATE(send_buffer,recv_buffer)

  END SUBROUTINE gather_i_2d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE allgather_r_1d_deblock(in_array, out_array, fill_value, &
    &                               allgather_pattern)
    ! dimension (nproma, nblk)
    REAL(dp), INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    REAL(dp), INTENT(INOUT) :: out_array(:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_allgather_pattern), INTENT(IN) :: allgather_pattern

    REAL(dp), ALLOCATABLE :: send_buffer(:,:)
    REAL(dp), POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk, n_procs, comm
    INTEGER, ALLOCATABLE :: collector_buffer_sizes(:)

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(allgather_pattern%gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(allgather_pattern%gather_pattern%loc_index(:))
      idx = idx_no(allgather_pattern%gather_pattern%loc_index(i))
      blk = blk_no(allgather_pattern%gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_r=send_buffer, fill_value=fill_value,&
      gather_pattern=allgather_pattern%gather_pattern, &
      collector_buffer_r=collector_buffer)
    DEALLOCATE(send_buffer)
    IF (allgather_pattern%intercomm /= MPI_COMM_NULL) THEN
      n_procs = p_comm_remote_size(allgather_pattern%intercomm)
      comm = allgather_pattern%intercomm
    ELSE
      n_procs = p_n_work
      comm = p_comm_work
    END IF
    ALLOCATE(collector_buffer_sizes(n_procs))
    CALL p_allgather(SIZE(collector_buffer, 2), collector_buffer_sizes, &
         comm=comm)
    IF (SIZE(out_array, 1) < SUM(collector_buffer_sizes)) &
      CALL finish("allgather_r_1d_deblock", "invalid out_array size")
#if defined (__SX__) || defined (__NEC_VH__)
    ! Workaround for occasional segfaults
    CALL p_barrier(comm)
#endif
    CALL p_allgatherv(collector_buffer(1,:), out_array, collector_buffer_sizes,&
      &               comm=comm)

    DEALLOCATE(collector_buffer_sizes, collector_buffer)
  END SUBROUTINE allgather_r_1d_deblock


  !-------------------------------------------------------------------------


  SUBROUTINE allgather_i_1d_deblock(in_array, out_array, fill_value, &
    &                               allgather_pattern)
    ! dimension (nproma, nblk)
    INTEGER, INTENT(IN) :: in_array(:,:)
    ! dimension (global length); only required on root
    INTEGER, INTENT(INOUT) :: out_array(:)
    INTEGER, INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_allgather_pattern), INTENT(IN) :: allgather_pattern

    INTEGER, ALLOCATABLE :: send_buffer(:,:)
    INTEGER, POINTER :: collector_buffer(:,:)
    INTEGER :: i, num_send_points, idx, blk, n_procs, comm
    INTEGER, ALLOCATABLE :: collector_buffer_sizes(:)

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    num_send_points = SUM(allgather_pattern%gather_pattern%collector_send_size(:))
    ALLOCATE(send_buffer(1, num_send_points))

    DO i = 1, SIZE(allgather_pattern%gather_pattern%loc_index(:))
      idx = idx_no(allgather_pattern%gather_pattern%loc_index(i))
      blk = blk_no(allgather_pattern%gather_pattern%loc_index(i))
      send_buffer(1,i) = in_array(idx, blk)
    END DO

    CALL two_phase_gather_first(send_buffer_i=send_buffer, fill_value=fill_value,&
      gather_pattern=allgather_pattern%gather_pattern, &
      collector_buffer_i=collector_buffer)
    DEALLOCATE(send_buffer)
    IF (allgather_pattern%intercomm /= MPI_COMM_NULL) THEN
      n_procs = p_comm_remote_size(allgather_pattern%intercomm)
      comm = allgather_pattern%intercomm
    ELSE
      n_procs = p_n_work
      comm = p_comm_work
    END IF
    ALLOCATE(collector_buffer_sizes(n_procs))
    CALL p_allgather(SIZE(collector_buffer, 2), collector_buffer_sizes, &
         comm=comm)
    IF (SIZE(out_array, 1) < SUM(collector_buffer_sizes)) &
      CALL finish("allgather_i_1d_deblock", "invalid out_array size")
    CALL p_allgatherv(collector_buffer(1,:), out_array, collector_buffer_sizes,&
      &               comm=comm)

    DEALLOCATE(collector_buffer_sizes, collector_buffer)
  END SUBROUTINE allgather_i_1d_deblock


  !-------------------------------------------------------------------------

  SUBROUTINE two_phase_gather_first_param_setup(use_fill_value,  &
       num_send_per_process, send_displ,                         &
       num_recv_per_process, recv_displ,                         &
       collector_buffer_nofill_size, collector_buffer_fill_size, &
       gather_pattern, fill_value_is_present)
    LOGICAL, INTENT(out) :: use_fill_value
    INTEGER, INTENT(out) :: num_send_per_process(p_n_work), send_displ(p_n_work), &
         num_recv_per_process(p_n_work), recv_displ(p_n_work), &
         collector_buffer_nofill_size, collector_buffer_fill_size
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    LOGICAL, INTENT(in) :: fill_value_is_present
    INTEGER :: i, collector_idx, num_coll, num_points_per_coll, &
         accum_send, accum_recv

    num_coll = SIZE(gather_pattern%collector_pes)
    collector_idx = -1
    DO i = 1, num_coll
      IF (gather_pattern%collector_pes(i) == p_pe_work) THEN
        collector_idx = i
        EXIT
      END IF
    END DO

    IF (collector_idx /= -1) THEN
      IF (fill_value_is_present) THEN
        num_points_per_coll = (gather_pattern%global_size + num_coll - 1) &
          &                   / num_coll
        collector_buffer_fill_size = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (collector_idx- 1)))
        use_fill_value = (collector_buffer_fill_size /= &
          &               gather_pattern%collector_size(collector_idx)) .AND. &
          &              (collector_buffer_fill_size > 0)
      ELSE
        collector_buffer_fill_size = gather_pattern%collector_size(collector_idx)
        use_fill_value = .FALSE.
      END IF
      collector_buffer_nofill_size = &
        gather_pattern%collector_size(collector_idx)
    ELSE
      use_fill_value = .FALSE.
      collector_buffer_nofill_size = 0
      collector_buffer_fill_size = 0
    END IF

    num_send_per_process(:) = 0
    num_send_per_process(gather_pattern%collector_pes(:)+1) = &
      gather_pattern%collector_send_size(:)
    num_recv_per_process(:) = 0
    num_recv_per_process(gather_pattern%recv_pes(:)+1) = &
      gather_pattern%recv_size(:)
    accum_send = 0
    accum_recv = 0
    DO i = 1, p_n_work
      send_displ(i) = accum_send
      accum_send = accum_send + num_send_per_process(i)
      recv_displ(i) = accum_recv
      accum_recv = accum_recv + num_recv_per_process(i)
    END DO

  END SUBROUTINE two_phase_gather_first_param_setup

  SUBROUTINE two_phase_gather_first_r(send_buffer_r, fill_value,                &
    &                                 gather_pattern, collector_buffer_r)
    ! dimension (:, length), prepared according to gather pattern
    REAL(dp), INTENT(IN) :: send_buffer_r(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(dp), POINTER, INTENT(OUT) :: collector_buffer_r(:,:)

    REAL(dp), POINTER :: collector_buffer_nofill_r(:,:)
    REAL(dp), POINTER :: collector_buffer_fill_r(:,:)
    INTEGER :: num_send_per_process(p_n_work), send_displ(p_n_work)
    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: collector_buffer_nofill_size, collector_buffer_fill_size
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    CALL two_phase_gather_first_param_setup(use_fill_value, &
         num_send_per_process, send_displ, &
         num_recv_per_process, recv_displ, &
         collector_buffer_nofill_size, collector_buffer_fill_size, &
         gather_pattern, PRESENT(fill_value))

    ALLOCATE(collector_buffer_nofill_r(SIZE(send_buffer_r, 1), &
      &      collector_buffer_nofill_size))

    CALL p_alltoallv(send_buffer_r, num_send_per_process, send_displ, &
      &              collector_buffer_nofill_r, num_recv_per_process, &
      &              recv_displ, p_comm_work)

    ! reorder collector_buffer
    collector_buffer_nofill_r(:,:) = &
      collector_buffer_nofill_r(:, gather_pattern%recv_buffer_reorder)
    IF (use_fill_value) THEN
      ALLOCATE(collector_buffer_fill_r(SIZE(send_buffer_r, 1), &
        &      collector_buffer_fill_size))
      collector_buffer_fill_r = fill_value
      collector_buffer_fill_r(:,gather_pattern%recv_buffer_reorder_fill) = &
        collector_buffer_nofill_r
      DEALLOCATE(collector_buffer_nofill_r)
      collector_buffer_r => collector_buffer_fill_r
    ELSE
      collector_buffer_r => collector_buffer_nofill_r
    END IF
  END SUBROUTINE two_phase_gather_first_r


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_first_s(send_buffer_r, fill_value, &
    &                                 gather_pattern, collector_buffer_r)
    ! dimension (:, length), prepared according to gather pattern
    REAL(sp), INTENT(IN) :: send_buffer_r(:,:)
    REAL(sp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(sp), POINTER, INTENT(OUT) :: collector_buffer_r(:,:)

    REAL(sp), POINTER :: collector_buffer_nofill_r(:,:)
    REAL(sp), POINTER :: collector_buffer_fill_r(:,:)
    INTEGER :: num_send_per_process(p_n_work), send_displ(p_n_work)
    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: collector_buffer_nofill_size, collector_buffer_fill_size
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    CALL two_phase_gather_first_param_setup(use_fill_value, &
         num_send_per_process, send_displ, &
         num_recv_per_process, recv_displ, &
         collector_buffer_nofill_size, collector_buffer_fill_size, &
         gather_pattern, PRESENT(fill_value))

    ALLOCATE(collector_buffer_nofill_r(SIZE(send_buffer_r, 1), &
      &      collector_buffer_nofill_size))

    CALL p_alltoallv(send_buffer_r, num_send_per_process, send_displ, &
      &              collector_buffer_nofill_r, num_recv_per_process, &
      &              recv_displ, p_comm_work)

    ! reorder collector_buffer
    collector_buffer_nofill_r(:,:) = &
      collector_buffer_nofill_r(:, gather_pattern%recv_buffer_reorder)
    IF (use_fill_value) THEN
      ALLOCATE(collector_buffer_fill_r(SIZE(send_buffer_r, 1), &
        &      collector_buffer_fill_size))
      collector_buffer_fill_r = fill_value
      collector_buffer_fill_r(:,gather_pattern%recv_buffer_reorder_fill) = &
        collector_buffer_nofill_r
      DEALLOCATE(collector_buffer_nofill_r)
      collector_buffer_r => collector_buffer_fill_r
    ELSE
      collector_buffer_r => collector_buffer_nofill_r
    END IF
  END SUBROUTINE two_phase_gather_first_s


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_first_i(send_buffer_i, fill_value, &
    &                                 gather_pattern, collector_buffer_i)
    ! dimension (:, length), prepared according to gather pattern
    INTEGER, INTENT(IN) :: send_buffer_i(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    INTEGER, POINTER, INTENT(OUT) :: collector_buffer_i(:,:)

    INTEGER, POINTER :: collector_buffer_nofill_i(:,:)
    INTEGER, POINTER :: collector_buffer_fill_i(:,:)
    INTEGER :: num_send_per_process(p_n_work), send_displ(p_n_work)
    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    INTEGER :: collector_buffer_nofill_size, collector_buffer_fill_size
    LOGICAL :: use_fill_value

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    CALL two_phase_gather_first_param_setup(use_fill_value, &
         num_send_per_process, send_displ, &
         num_recv_per_process, recv_displ, &
         collector_buffer_nofill_size, collector_buffer_fill_size, &
         gather_pattern, PRESENT(fill_value))

    ALLOCATE(collector_buffer_nofill_i(SIZE(send_buffer_i, 1), &
      &      collector_buffer_nofill_size))

    CALL p_alltoallv(send_buffer_i, num_send_per_process, send_displ, &
      &              collector_buffer_nofill_i, num_recv_per_process, &
      &              recv_displ, p_comm_work)

    ! reorder collector_buffer
    collector_buffer_nofill_i(:,:) = &
      collector_buffer_nofill_i(:, gather_pattern%recv_buffer_reorder)
    IF (use_fill_value) THEN
      ALLOCATE(collector_buffer_fill_i(SIZE(send_buffer_i, 1), &
        &      collector_buffer_fill_size))
      collector_buffer_fill_i = fill_value
      collector_buffer_fill_i(:,gather_pattern%recv_buffer_reorder_fill) = &
        collector_buffer_nofill_i
      DEALLOCATE(collector_buffer_nofill_i)
      collector_buffer_i => collector_buffer_fill_i
    ELSE
      collector_buffer_i => collector_buffer_nofill_i
    END IF
  END SUBROUTINE two_phase_gather_first_i


  !-------------------------------------------------------------------------

  SUBROUTINE two_phase_gather_second_param_setup(&
       num_recv_per_process, recv_displ,         &
       gather_pattern, fill_value_is_present)
    INTEGER, INTENT(out) :: num_recv_per_process(p_n_work), recv_displ(p_n_work)
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    LOGICAL, INTENT(in) :: fill_value_is_present

    INTEGER :: i, num_coll, collector_idx, num_points_per_coll, &
         collector_buffer_fill_size, recv_accum
    LOGICAL :: use_fill_value

    num_coll = SIZE(gather_pattern%collector_pes)
    collector_idx = -1
    DO i = 1, num_coll
      IF (gather_pattern%collector_pes(i) == p_pe_work) THEN
        collector_idx = i
        EXIT
      END IF
    END DO

    IF (collector_idx /= -1 .AND. fill_value_is_present) THEN
      num_points_per_coll = (gather_pattern%global_size + num_coll - 1) &
        &                   / num_coll
      collector_buffer_fill_size = &
        & MIN(num_points_per_coll, &
        &     gather_pattern%global_size - &
        &       MAX(0,num_points_per_coll * (collector_idx - 1)))
      use_fill_value = (collector_buffer_fill_size /= &
           &              gather_pattern%collector_size(collector_idx)) &
           &        .AND. collector_buffer_fill_size > 0
    ELSE
      use_fill_value = .FALSE.
    END IF

    num_recv_per_process(:) = 0
    IF (p_pe_work == process_mpi_root_id) THEN
      IF (use_fill_value) THEN
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          MIN(num_points_per_coll, &
          & gather_pattern%global_size - &
          & MAX(0,num_points_per_coll * (/(i, i = 0, num_coll - 1)/)))
      ELSE
        num_recv_per_process(gather_pattern%collector_pes(:)+1) = &
          gather_pattern%collector_size(:)
      END IF
    END IF
    recv_accum = 0
    DO i = 1, p_n_work
      recv_displ(i) = recv_accum
      recv_accum = recv_accum + num_recv_per_process(i)
    END DO
  END SUBROUTINE two_phase_gather_second_param_setup

  SUBROUTINE two_phase_gather_second_r(recv_buffer_r, fill_value, &
    &                                  gather_pattern, collector_buffer_r)
    ! dimension (:, global length); only required on root
    REAL(dp), INTENT(INOUT) :: recv_buffer_r(:,:)
    REAL(dp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(dp), POINTER, INTENT(INOUT) :: collector_buffer_r(:,:)

    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    CALL two_phase_gather_second_param_setup(&
         num_recv_per_process, recv_displ,   &
         gather_pattern, PRESENT(fill_value))

    CALL p_gatherv(collector_buffer_r, SIZE(collector_buffer_r, 2), &
      &            recv_buffer_r, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)
    DEALLOCATE(collector_buffer_r)
  END SUBROUTINE two_phase_gather_second_r


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_second_s(recv_buffer_r, fill_value,                &
    &                                  gather_pattern, collector_buffer_r)
    ! dimension (:, global length); only required on root
    REAL(sp), INTENT(INOUT) :: recv_buffer_r(:,:)
    REAL(sp), INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    REAL(sp), POINTER, INTENT(INOUT) :: collector_buffer_r(:,:)

    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    CALL two_phase_gather_second_param_setup(&
         num_recv_per_process, recv_displ,   &
         gather_pattern, PRESENT(fill_value))

    CALL p_gatherv(collector_buffer_r, SIZE(collector_buffer_r, 2), &
      &            recv_buffer_r, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)
    DEALLOCATE(collector_buffer_r)
  END SUBROUTINE two_phase_gather_second_s


  !-------------------------------------------------------------------------


  SUBROUTINE two_phase_gather_second_i(recv_buffer_i, fill_value,       &
    &                                  gather_pattern,                  &
    &                                  collector_buffer_i)
    ! dimension (:, global length); only required on root
    INTEGER, INTENT(INOUT) :: recv_buffer_i(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: fill_value ! if provided missing values will
    ! be replaced with this value
    ! if not provided all valid
    ! points will be packed to the
    ! front of the array
    TYPE(t_comm_gather_pattern), INTENT(IN) :: gather_pattern
    INTEGER, POINTER, INTENT(INOUT) :: collector_buffer_i(:,:)

    INTEGER :: num_recv_per_process(p_n_work), recv_displ(p_n_work)

    !
    ! OPENACC:  GPU execution assumes that all information is now on the host
    !

    CALL two_phase_gather_second_param_setup(&
         num_recv_per_process, recv_displ,   &
         gather_pattern, PRESENT(fill_value))

    CALL p_gatherv(collector_buffer_i, SIZE(collector_buffer_i, 2), &
      &            recv_buffer_i, num_recv_per_process, recv_displ, &
      &            process_mpi_root_id, p_comm_work)
    DEALLOCATE(collector_buffer_i)
  END SUBROUTINE two_phase_gather_second_i


  !-----------------------------------------------------------------------------
  !> Factory method for t_ScatterPattern. Destroy with deleteScatterPattern().
  !-----------------------------------------------------------------------------
  FUNCTION makeScatterPattern(jg, loc_arr_len, glb_index, communicator)
    USE mo_scatter_pattern_scatter
    IMPLICIT NONE
    CLASS(t_ScatterPattern), POINTER :: makeScatterPattern
    INTEGER, VALUE :: jg, loc_arr_len, communicator
    INTEGER, INTENT(IN) :: glb_index(:)

    CHARACTER(*), PARAMETER :: routine = modname//":makeScatterPattern"
    INTEGER :: ierr

    ALLOCATE(t_ScatterPatternScatter::makeScatterPattern, stat = ierr)
    IF(ierr /= SUCCESS) CALL finish(routine, "error allocating memory")
    CALL makeScatterPattern%construct(jg, loc_arr_len, glb_index, communicator)
  END FUNCTION makeScatterPattern

END MODULE mo_communication


