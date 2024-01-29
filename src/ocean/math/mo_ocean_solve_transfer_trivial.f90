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

#if (defined(_OPENMP) && defined(OCE_SOLVE_OMP))
#include "omp_definitions.inc"
#endif

MODULE mo_ocean_solve_trivial_transfer
  USE mo_kind, ONLY: wp, sp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_aux, ONLY: solve_cell, solve_edge, solve_vert
  USE mo_model_domain, ONLY: t_patch
  USE mo_mpi, ONLY: p_pe_work, p_comm_work, my_process_is_mpi_parallel
  USE mo_parallel_config, ONLY: nproma
  USE mo_timer, ONLY: timer_start, timer_stop, new_timer
  USE mo_run_config, ONLY: ltimer
  USE mo_communication, ONLY: t_comm_pattern, exchange_data
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

! provides extended communication / transfer infrastructure object derived from abstract t_transfer - type
! to be used by solvers
! trivial transfer : group of solver-PEs is same as group od solver-PEs
! arrays are just locally copied... (and converted between different real-kinds, if necessary)

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_trivial_transfer

  TYPE, EXTENDS(t_transfer) :: t_trivial_transfer
    PRIVATE
    CLASS(t_comm_pattern), POINTER :: comm_pat_sync
  CONTAINS
! overrides for deferred interfaces from parenting abstract type t_transfer
    PROCEDURE, PUBLIC :: construct => trivial_transfer_construct
    PROCEDURE, PUBLIC :: destruct => trivial_transfer_destruct
    PROCEDURE :: into_2d_wp => trivial_transfer_into_2d_wp
    PROCEDURE :: into_2d_wp_2 => trivial_transfer_into_2d_wp_2
    PROCEDURE :: into_3d_wp => trivial_transfer_into_3d_wp
    PROCEDURE :: into_idx => trivial_transfer_into_idx
    PROCEDURE :: into_once_2d_wp => trivial_transfer_into_once_2d_wp
    PROCEDURE :: into_once_3d_wp => trivial_transfer_into_once_3d_wp
    PROCEDURE :: into_once_idx => trivial_transfer_into_once_idx
    PROCEDURE :: out_2d_wp => trivial_transfer_out_2d_wp
    PROCEDURE :: bcst_1d_wp => trivial_transfer_bcst_1d_wp
    PROCEDURE :: bcst_1d_i => trivial_transfer_bcst_1d_i
    PROCEDURE :: sync_2d_wp => trivial_transfer_sync_2d_wp
    PROCEDURE :: sync_2d_sp => trivial_transfer_sync_2d_sp
  END TYPE t_trivial_transfer

  CHARACTER(LEN=*), PARAMETER :: module_name = "mo_ocean_solve_trivial_transfer"

CONTAINS

  SUBROUTINE trivial_transfer_construct(this, sync_type, patch_2d)
    CLASS(t_trivial_transfer), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: sync_type
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
      & "::trivial_transfer_construct()"
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: gID_tmp
    INTEGER :: iidx, iblk, nidx
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: gID_tmp
#endif

    IF (this%is_init) RETURN
    IF (ltimer) THEN
      this%timer_init = new_timer("triv-t init")
      CALL timer_start(this%timer_init)
    END IF
    this%is_solver_pe = .true.
    this%is_leader_pe = (0 .EQ. p_pe_work)
    this%comm = p_comm_work
    SELECT CASE(sync_type)
    CASE(solve_cell)
      this%nblk_a = patch_2D%alloc_cell_blocks
      this%nblk = patch_2d%cells%in_domain%end_block
      this%nidx = nproma
      this%nidx_l = nproma
      this%nidx_e = patch_2d%cells%in_domain%end_index
      this%glb_idx_loc => patch_2d%cells%decomp_info%glb_index
      this%comm_pat_sync => patch_2d%comm_pat_c
    CASE(solve_edge)
      this%nblk_a = SIZE(patch_2D%edges%decomp_info%owner_mask, 2)
      this%nblk = patch_2d%edges%in_domain%end_block
      this%nidx = nproma
      this%nidx_l = nproma
      this%nidx_e = patch_2d%edges%in_domain%end_index
      this%glb_idx_loc => patch_2d%edges%decomp_info%glb_index
      this%comm_pat_sync => patch_2d%comm_pat_e
    CASE(solve_vert)
      this%nblk_a = SIZE(patch_2D%verts%decomp_info%owner_mask, 2)
      this%nblk = patch_2d%verts%in_domain%end_block
      this%nidx = nproma
      this%nidx_l = nproma
      this%nidx_e = patch_2d%verts%in_domain%end_index
      this%glb_idx_loc => patch_2d%verts%decomp_info%glb_index
      this%comm_pat_sync => patch_2d%comm_pat_v
    CASE DEFAULT
      CALL finish(routine, "syncing scheme not recognized")
    END SELECT
    ALLOCATE(gID_tmp(this%nidx_l,this%nblk_a))
    DO iblk = 1, this%nblk
      nidx = MERGE(this%nidx_l, this%nidx_e, this%nblk .NE. iblk)
      DO iidx = 1, nidx
        gID_tmp(iidx, iblk) = this%globalID_loc(iidx, iblk)
      END DO
      IF (nidx .NE. this%nidx_l) &
        & gID_tmp(nidx+1:this%nidx_l, iblk) = -1
    END DO
    DO iblk = this%nblk + 1, this%nblk_a
      gID_tmp(1:this%nidx_l, iblk) = -1
    ENDDO
    IF(my_process_is_mpi_parallel()) &
      & CALL exchange_data(this%comm_pat_sync, gID_tmp)
    NULLIFY(this%glb_idx_loc)
    ALLOCATE(this%glb_idx_loc(this%nidx_l * this%nblk_a))
    DO iblk = 1, this%nblk_a
      DO iidx = 1, this%nidx_l
        this%glb_idx_loc(iidx + (iblk - 1) * this%nidx_l) = gID_tmp(iidx, iblk)
      END DO
    END DO
    DEALLOCATE(gID_tmp)
    this%ngid_a_l = SIZE(this%glb_idx_loc)
    this%glb_idx_cal => this%glb_idx_loc
    IF (ltimer) THEN
      this%timer_glob_sum = new_timer("triv-t glb_sum")
      this%timer_sync = new_timer("triv-t sync")
      this%timer_in(1) = new_timer("triv-t in(solve)")
      this%timer_in(2) = new_timer("triv-t in(lhs)")
      this%timer_in(3) = new_timer("triv-t in(init solve)")
      this%timer_in(4) = new_timer("triv-t in(init lhs)")
      this%timer_out = new_timer("triv-t out")
      CALL timer_stop(this%timer_init)
    END IF
    this%is_init = .true.
  END SUBROUTINE trivial_transfer_construct

  SUBROUTINE trivial_transfer_destruct(this)
    CLASS(t_trivial_transfer), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%glb_idx_loc)) DEALLOCATE(this%glb_idx_loc)
    NULLIFY(this%glb_idx_loc, this%glb_idx_cal)
    this%is_init = .false.
  END SUBROUTINE trivial_transfer_destruct

  SUBROUTINE trivial_transfer_into_once_2d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(data_in) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_in(tt))
    IF (.NOT.ALLOCATED(data_out)) THEN
        ALLOCATE(data_out(SIZE(data_in, 1), SIZE(data_in, 2)))
        !$ACC ENTER DATA CREATE(data_out) IF(lzacc)
    END IF
!ICON_OMP PARALLEL WORKSHARE
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    data_out(:,:) = data_in(:,:)
    !$ACC END KERNELS
    !$ACC WAIT(1)
!ICON_OMP END PARALLEL WORKSHARE
    IF (ltimer) CALL timer_stop(this%timer_in(tt))

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_into_once_2d_wp

  SUBROUTINE trivial_transfer_into_once_3d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (.NOT.ALLOCATED(data_out)) &
      & ALLOCATE(data_out(SIZE(data_in, 1), SIZE(data_in, 2), SIZE(data_in, 3)))

    CALL this%into(data_in, data_out, tt, lacc=lzacc)

  END SUBROUTINE trivial_transfer_into_once_3d_wp

  SUBROUTINE trivial_transfer_into_once_idx(this, data_in_idx, data_in_blk, &
     &  data_out_idx, data_out_blk, tt, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in_idx, data_in_blk
    INTEGER, INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: &
      & data_out_idx, data_out_blk
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (.NOT.ALLOCATED(data_out_idx)) &
      & ALLOCATE(data_out_idx(SIZE(data_in_idx, 1), &
        & SIZE(data_in_idx, 2), SIZE(data_in_idx, 3)), &
        & data_out_blk(SIZE(data_in_blk, 1), &
        & SIZE(data_in_blk, 2), SIZE(data_in_blk, 3)))

    CALL this%into(data_in_idx, data_in_blk, &
      & data_out_idx, data_out_blk, tt, lacc=lzacc)

  END SUBROUTINE trivial_transfer_into_once_idx

  SUBROUTINE trivial_transfer_into_2d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(data_in) &
    !$ACC   COPYOUT(data_out) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_in(tt))
!ICON_OMP PARALLEL WORKSHARE
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    data_out(:,:) = data_in(:,:)
    !$ACC END KERNELS
    !$ACC WAIT(1)
!ICON_OMP END PARALLEL WORKSHARE
    IF (ltimer) CALL timer_stop(this%timer_in(tt))

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_into_2d_wp

  SUBROUTINE trivial_transfer_into_2d_wp_2(this, di1, do1, di2, do2, tt, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: di1, di2
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: do1, do2
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(di1, di2) &
    !$ACC   COPY(do1, do2) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_in(tt))
!ICON_OMP PARALLEL WORKSHARE
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    do1(:,:) = di1(:,:)
    do2(:,:) = di2(:,:)
    !$ACC END KERNELS
    !$ACC WAIT(1)
!ICON_OMP END PARALLEL WORKSHARE
    IF (ltimer) CALL timer_stop(this%timer_in(tt))

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_into_2d_wp_2

  SUBROUTINE trivial_transfer_into_3d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:,:), CONTIGUOUS :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    INTEGER :: i

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(data_in) &
    !$ACC   COPYOUT(data_out) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_in(tt))
#ifdef _OPENMP
!$OMP PARALLEL DO SCHEDULE(STATIC)
#endif
    DO i = 1, SIZE(data_in, 3)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      data_out(:,:,i) = data_in(:,:,i)
      !$ACC END KERNELS
      !$ACC WAIT(1)
    END DO
#ifdef _OPENMP
!$OMP END PARALLEL DO
#endif
    IF (ltimer) CALL timer_stop(this%timer_in(tt))

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_into_3d_wp

  SUBROUTINE trivial_transfer_into_idx(this, data_in_idx, data_in_blk, &
     &  data_out_idx, data_out_blk, tt, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in_blk, data_in_idx
    INTEGER, INTENT(OUT), DIMENSION(:,:,:), CONTIGUOUS :: data_out_blk, data_out_idx
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    INTEGER :: i
    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(data_in_idx, data_in_blk) &
    !$ACC   COPYOUT(data_out_idx, data_out_blk) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_in(tt))
#ifdef _OPENMP
!$OMP PARALLEL DO SCHEDULE(STATIC)
#endif
    DO i = 1, SIZE(data_in_idx, 3)
      !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
      data_out_idx(:,:,i) = data_in_idx(:,:,i)
      data_out_blk(:,:,i) = data_in_blk(:,:,i)
      !$ACC END KERNELS
      !$ACC WAIT(1)
    END DO
#ifdef _OPENMP
!$OMP END PARALLEL DO
#endif
    IF (ltimer) CALL timer_stop(this%timer_in(tt))

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_into_idx

  SUBROUTINE trivial_transfer_out_2d_wp(this, data_in, data_out, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: data_out
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(data_in) &
    !$ACC   COPY(data_out) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_out)
!ICON_OMP PARALLEL WORKSHARE
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    data_out(:,:) = data_in(:,:)
    !$ACC END KERNELS
    !$ACC WAIT(1)
!ICON_OMP END PARALLEL WORKSHARE
    IF (ltimer) CALL timer_stop(this%timer_out)

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_out_2d_wp

  SUBROUTINE trivial_transfer_bcst_1d_wp(this, data_in, data_out, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:), CONTIGUOUS :: data_out
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(data_in) &
    !$ACC   COPYOUT(data_out) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_out)
!ICON_OMP PARALLEL WORKSHARE
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    data_out(:) = data_in(:)
    !$ACC END KERNELS
    !$ACC WAIT(1)
!ICON_OMP END PARALLEL WORKSHARE
    IF (ltimer) CALL timer_stop(this%timer_out)

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_bcst_1d_wp

  SUBROUTINE trivial_transfer_bcst_1d_i(this, data_in, data_out, lacc)
    CLASS(t_trivial_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN), DIMENSION(:), CONTIGUOUS :: data_in
    INTEGER, INTENT(OUT), DIMENSION(:), CONTIGUOUS :: data_out
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    !$ACC DATA COPYIN(data_in) &
    !$ACC   COPYOUT(data_out) IF(lzacc)

    IF (ltimer) CALL timer_start(this%timer_out)
!ICON_OMP PARALLEL WORKSHARE
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
    data_out(:) = data_in(:)
    !$ACC END KERNELS
    !$ACC WAIT(1)
!ICON_OMP END PARALLEL WORKSHARE
    IF (ltimer) CALL timer_stop(this%timer_out)

    !$ACC END DATA

  END SUBROUTINE trivial_transfer_bcst_1d_i

  SUBROUTINE trivial_transfer_sync_2d_wp(this, data_inout)
    CLASS(t_trivial_transfer), INTENT(INOUT) :: this
    REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout

    IF (ltimer) CALL timer_start(this%timer_sync)
    IF(my_process_is_mpi_parallel()) &
      CALL exchange_data(this%comm_pat_sync, data_inout)
    IF (ltimer) CALL timer_stop(this%timer_sync)
  END SUBROUTINE trivial_transfer_sync_2d_wp

  SUBROUTINE trivial_transfer_sync_2d_sp(this, data_inout)
    CLASS(t_trivial_transfer), INTENT(INOUT) :: this
    REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout

    IF (ltimer) CALL timer_start(this%timer_sync)
    IF(my_process_is_mpi_parallel()) &
      CALL exchange_data(this%comm_pat_sync, data_inout)
    IF (ltimer) CALL timer_stop(this%timer_sync)
  END SUBROUTINE trivial_transfer_sync_2d_sp

END MODULE mo_ocean_solve_trivial_transfer
