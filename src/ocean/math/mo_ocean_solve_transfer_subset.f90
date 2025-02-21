! provides extended communication / transfer infrastructure object
! derived from abstract t_transfer - type to be used by solvers
!
! trivial transfer : group of solver-PEs is same as group od
! solver-PEs arrays are just locally copied... (and converted between
! different real-kinds, if necessary)
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





MODULE mo_ocean_solve_subset_transfer
  USE mo_kind, ONLY: wp, sp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_aux, ONLY: solve_trans_scatter, &
    & solve_trans_compact, solve_cell, solve_edge, solve_vert
  USE mo_model_domain, ONLY: t_patch
  USE mo_mpi, ONLY: p_n_work, p_pe_work, p_comm_work, p_sum, p_int, &
    & p_bcast, my_process_is_mpi_parallel
  USE mo_parallel_config, ONLY: nproma
  USE mo_timer, ONLY: timer_start, timer_stop, new_timer
  USE mo_run_config, ONLY: ltimer
  USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, &
    & init_glb2loc_index_lookup, set_inner_glb_index, deallocate_glb2loc_index_lookup
  USE mo_communication, ONLY: t_comm_pattern, delete_comm_pattern, &
    & exchange_data, exchange_data_mult
  USE mo_communication_factory, ONLY: setup_comm_pattern
  USE mo_fortran_tools, ONLY: set_acc_host_or_device

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_subset_transfer

  TYPE, EXTENDS(t_transfer) :: t_subset_transfer
    PRIVATE
    INTEGER :: solve_buddy, nblk_l, nblk_a_l, nidx_e_l, ngid_a
    CLASS(t_comm_pattern), POINTER :: cpat_in => NULL(), &
      & cpat_in2 => NULL(), cpat_out => NULL(), cpat_sync => NULL()
  CONTAINS
! overrides for deferred interfaces from parenting abstract type t_transfer
    PROCEDURE, PUBLIC :: construct => subset_transfer_construct
    PROCEDURE, PUBLIC :: destruct => subset_transfer_destruct
    PROCEDURE :: into_2d_wp => subset_transfer_into_2d_wp
    PROCEDURE :: into_2d_wp_2 => subset_transfer_into_2d_wp_2
    PROCEDURE :: into_3d_wp => subset_transfer_into_3d_wp
    PROCEDURE :: into_idx => subset_transfer_into_idx
    PROCEDURE :: into_once_2d_wp => subset_transfer_into_once_2d_wp
    PROCEDURE :: into_once_3d_wp => subset_transfer_into_once_3d_wp
    PROCEDURE :: into_once_idx => subset_transfer_into_once_idx
    PROCEDURE :: out_2d_wp => subset_transfer_out_2d_wp
    PROCEDURE :: bcst_1d_wp => subset_transfer_bcst_1d_wp
    PROCEDURE :: bcst_1d_i => subset_transfer_bcst_1d_i
    PROCEDURE :: sync_2d_wp => subset_transfer_sync_2d_wp
    PROCEDURE :: sync_2d_sp => subset_transfer_sync_2d_sp
  END TYPE t_subset_transfer

  CHARACTER(LEN=*), PARAMETER :: module_name = "mo_ocean_solve_subset_transfer"

CONTAINS

  SUBROUTINE subset_transfer_construct(this, sync_type, patch_2d, redfac, mode)
    CLASS(t_subset_transfer), INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: sync_type, redfac, mode
    TYPE(t_patch), POINTER :: patch_2d
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
      & "::subset_transfer_construct()"

    CALL finish(routine, "subset transfer does not make sense without MPI")
  END SUBROUTINE subset_transfer_construct

  SUBROUTINE subset_transfer_destruct(this)
    CLASS(t_subset_transfer), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%cpat_in)) THEN
      CALL delete_comm_pattern(this%cpat_in)
      CALL delete_comm_pattern(this%cpat_in2)
      CALL delete_comm_pattern(this%cpat_out)
    END IF
    IF (ASSOCIATED(this%cpat_sync)) &
      & CALL delete_comm_pattern(this%cpat_sync)
    IF (ASSOCIATED(this%glb_idx_loc)) DEALLOCATE(this%glb_idx_loc)
    IF (ASSOCIATED(this%glb_idx_cal)) DEALLOCATE(this%glb_idx_cal)
    NULLIFY(this%glb_idx_loc, this%glb_idx_cal)
    this%is_init = .false.
  END SUBROUTINE subset_transfer_destruct

  SUBROUTINE subset_transfer_into_once_2d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_into_once_2d_wp()"

    CALL set_acc_host_or_device(lzacc, lacc)


    IF (.NOT.ALLOCATED(data_out)) THEN
      ALLOCATE(data_out(this%nidx, this%nblk_a))
      IF (this%is_solver_pe) data_out(:, this%nblk_a) = 0._wp
    END IF
    CALL this%into(data_in, data_out, tt)
  END SUBROUTINE subset_transfer_into_once_2d_wp

  SUBROUTINE subset_transfer_into_once_3d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_into_once_3d_wp()"

    CALL set_acc_host_or_device(lzacc, lacc)


    IF (.NOT.ALLOCATED(data_out)) THEN
      ALLOCATE(data_out(this%nidx, this%nblk, SIZE(data_in, 3)))
      IF (this%is_solver_pe) data_out(:, this%nblk, :) = 0._wp
    END IF
    CALL this%into(data_in, data_out, tt)
  END SUBROUTINE subset_transfer_into_once_3d_wp

  SUBROUTINE subset_transfer_into_once_idx(this, data_in_idx, data_in_blk, &
     &  data_out_idx, data_out_blk, tt, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in_idx, data_in_blk
    INTEGER, INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: &
      & data_out_idx, data_out_blk
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_into_once_idx()"

    CALL set_acc_host_or_device(lzacc, lacc)


    IF (.NOT.ALLOCATED(data_out_idx)) &
      & ALLOCATE(data_out_idx(this%nidx, this%nblk, SIZE(data_in_idx, 3)), &
        & data_out_blk(this%nidx, this%nblk, SIZE(data_in_blk, 3)))
    CALL this%into(data_in_idx, data_in_blk, data_out_idx, data_out_blk, tt)
  END SUBROUTINE subset_transfer_into_once_idx

  SUBROUTINE subset_transfer_into_2d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_into_2d_wp()"

    CALL set_acc_host_or_device(lzacc, lacc)


    IF (ltimer) CALL timer_start(this%timer_in(tt))
    CALL exchange_data(this%cpat_in, data_out, data_in)
    IF (ltimer) CALL timer_stop(this%timer_in(tt))
  END SUBROUTINE subset_transfer_into_2d_wp

  SUBROUTINE subset_transfer_into_2d_wp_2(this, di1, do1, di2, do2, tt, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: di1, di2
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: do1, do2
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: to, ti
    INTEGER :: i
    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_into_2d_wp_2()"


    CALL set_acc_host_or_device(lzacc, lacc)


    IF (ltimer) CALL timer_start(this%timer_in(tt))
    ALLOCATE(to(SIZE(do1, 1), 1, SIZE(do1, 2), 2), &
      & ti(SIZE(di1, 1), 1, SIZE(di1, 2), 2))
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO i = 1, SIZE(di1, 2)
      ti(:,1,i,1) = di1(:,i)
      ti(:,1,i,2) = di2(:,i)
    END DO
!ICON_OMP END PARALLEL DO
    IF (this%is_solver_pe) to(:, 1, SIZE(do1, 2), :) = 0._wp
    CALL exchange_data_mult(this%cpat_in, 2, 2, recv4d=to, send4d=ti)
    IF (this%is_solver_pe) THEN
      DO i = 1, 2
        IF (i .EQ. 1) THEN
          do1(:,:) = to(:,1,:,i)
        ELSE
          do2(:,:) = to(:,1,:,i)
        END IF
      END DO
    END IF
    DEALLOCATE(ti, to)
    IF (ltimer) CALL timer_stop(this%timer_in(tt))
  END SUBROUTINE subset_transfer_into_2d_wp_2

  SUBROUTINE subset_transfer_into_3d_wp(this, data_in, data_out, tt, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:,:), CONTIGUOUS :: data_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    INTEGER :: i, j, n3
    REAL(KIND=wp), DIMENSION(:,:,:,:), ALLOCATABLE :: to
    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_into_3d_wp()"


    CALL set_acc_host_or_device(lzacc, lacc)


    IF (ltimer) CALL timer_start(this%timer_in(tt))
    n3 = SIZE(data_in, 3)
    ALLOCATE(to(this%nidx, 1, this%nblk, n3))
    IF (this%is_solver_pe) to(:, 1, this%nblk, :) = 0._wp
    CALL exchange_data_mult(this%cpat_in2, n3, n3, recv4d=to, send4d= &
      & RESHAPE(data_in, (/SIZE(data_in, 1), 1, SIZE(data_in, 2), n3/)))
    IF (this%is_solver_pe) THEN
      DO i = 1, n3
        data_out(:,:,i) = to(:,1,:,i)
      END DO
    END IF
    DEALLOCATE(to)
    IF (ltimer) CALL timer_stop(this%timer_in(tt))
  END SUBROUTINE subset_transfer_into_3d_wp

  SUBROUTINE subset_transfer_into_idx(this, data_in_idx, data_in_blk, &
     &  data_out_idx, data_out_blk, tt, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in_blk, data_in_idx
    INTEGER, INTENT(OUT), DIMENSION(:,:,:), CONTIGUOUS :: data_out_blk, data_out_idx
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: glb_in, glb_out
    INTEGER, INTENT(IN) :: tt
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    INTEGER :: i, iblk, iidx, jblk, jidx, gid
    LOGICAL :: found, notfound
    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_into_idx()"


    CALL set_acc_host_or_device(lzacc, lacc)


    IF (ltimer) CALL timer_start(this%timer_in(tt))
    DO i = 1, SIZE(data_in_idx, 3)
      ALLOCATE(glb_in(this%nidx_l, this%nblk_l), &
        & glb_out(this%nidx, this%nblk))
!ICON_OMP PARALLEL DO PRIVATE(iidx) SCHEDULE(STATIC)
      DO iblk = 1, this%nblk_l
        DO iidx = 1, MERGE(this%nidx_l, this%nidx_e_l,iblk .LT. this%nblk_l)
          glb_in(iidx, iblk) = this%globalID_loc(data_in_idx(iidx, iblk, i), data_in_blk(iidx, iblk, i))
        END DO
      END DO
!ICON_OMP END PARALLEL DO
      glb_in(this%nidx_e_l + 1:, this%nblk_l) = -1
      CALL exchange_data(this%cpat_in2, glb_out, glb_in)
      IF (this%is_solver_pe) THEN
        DO iblk = 1, this%nblk
          DO iidx = 1, MERGE(this%nidx, this%nidx_e,iblk .LT. this%nblk)
            gid = glb_out(iidx, iblk)
            found = .FALSE.
            DO jblk = 1, this%nblk_a
              DO jidx = 1, this%nidx
                IF (this%globalID_cal(jidx, jblk) .EQ. gid) THEN
                  found = .TRUE.
                  data_out_idx(iidx, iblk, i) = jidx
                  data_out_blk(iidx, iblk, i) = jblk
                  EXIT
                END IF
              END DO
              IF (found) EXIT
            END DO
          END DO
        END DO
        data_out_idx(this%nidx_e + 1:, this%nblk, i) = data_out_idx(this%nidx_e, this%nblk, i)
        data_out_blk(this%nidx_e + 1:, this%nblk, i) = data_out_blk(this%nidx_e, this%nblk, i)
      END IF
      DEALLOCATE(glb_in, glb_out)
    END DO
    IF (ltimer) CALL timer_stop(this%timer_in(tt))
  END SUBROUTINE subset_transfer_into_idx

  SUBROUTINE subset_transfer_out_2d_wp(this, data_in, data_out, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:,:), CONTIGUOUS :: data_out
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_out_2d_wp()"

    CALL set_acc_host_or_device(lzacc, lacc)


    IF (ltimer) CALL timer_start(this%timer_out)
    CALL exchange_data(this%cpat_out, data_out, data_in)
    IF (ltimer) CALL timer_stop(this%timer_out)
  END SUBROUTINE subset_transfer_out_2d_wp

  SUBROUTINE subset_transfer_bcst_1d_wp(this, data_in, data_out, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    REAL(KIND=wp), INTENT(IN), DIMENSION(:), CONTIGUOUS :: data_in
    REAL(KIND=wp), INTENT(OUT), DIMENSION(:), CONTIGUOUS :: data_out
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_bcst_1d_wp()"

    CALL set_acc_host_or_device(lzacc, lacc)


    IF (ltimer) CALL timer_start(this%timer_out)
    IF (this%is_leader_pe) data_out(:) = data_in(:)
    CALL p_bcast(data_out, 0)
    IF (ltimer) CALL timer_stop(this%timer_out)
  END SUBROUTINE subset_transfer_bcst_1d_wp

  SUBROUTINE subset_transfer_bcst_1d_i(this, data_in, data_out, lacc)
    CLASS(t_subset_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN), DIMENSION(:), CONTIGUOUS :: data_in
    INTEGER, INTENT(OUT), DIMENSION(:), CONTIGUOUS :: data_out
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    LOGICAL :: lzacc
    CHARACTER(LEN=*), PARAMETER :: routine = module_name// &
                            & "::subset_transfer_bcst_1d_i()"

    CALL set_acc_host_or_device(lzacc, lacc)


    IF (ltimer) CALL timer_start(this%timer_out)
    IF (this%is_leader_pe) data_out(:) = data_in(:)
    CALL p_bcast(data_out, 0)
    IF (ltimer) CALL timer_stop(this%timer_out)
  END SUBROUTINE subset_transfer_bcst_1d_i

  SUBROUTINE subset_transfer_sync_2d_wp(this, data_inout)
    CLASS(t_subset_transfer), INTENT(INOUT) :: this
    REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout

    IF (ltimer) CALL timer_start(this%timer_sync)
    IF (my_process_is_mpi_parallel()) &
      & CALL exchange_data(this%cpat_sync, data_inout)
    IF (ltimer) CALL timer_stop(this%timer_sync)
  END SUBROUTINE subset_transfer_sync_2d_wp

  SUBROUTINE subset_transfer_sync_2d_sp(this, data_inout)
    CLASS(t_subset_transfer), INTENT(INOUT) :: this
    REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout

    IF (ltimer) CALL timer_start(this%timer_sync)
    IF (my_process_is_mpi_parallel()) &
      & CALL exchange_data(this%cpat_sync, data_inout)
    IF (ltimer) CALL timer_stop(this%timer_sync)
  END SUBROUTINE subset_transfer_sync_2d_sp

END MODULE mo_ocean_solve_subset_transfer
