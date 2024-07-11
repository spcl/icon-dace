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

#if (defined(_OPENMP) && defined(OCE_SOLVE_OMP))
#include "omp_definitions.inc"
#endif

MODULE mo_ocean_solve_subset_transfer
  USE mo_kind, ONLY: wp, sp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_aux, ONLY: solve_trans_scatter, &
    & solve_trans_compact, solve_cell, solve_edge, solve_vert
  USE mo_model_domain, ONLY: t_patch
  USE mo_mpi, ONLY: p_n_work, p_pe_work, p_comm_work, p_sum, p_int, &
    & p_bcast, my_process_is_mpi_parallel, i_am_accel_node
  USE mo_parallel_config, ONLY: nproma
  USE mo_timer, ONLY: timer_start, timer_stop, new_timer
  USE mo_run_config, ONLY: ltimer
  USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, &
    & init_glb2loc_index_lookup, set_inner_glb_index, deallocate_glb2loc_index_lookup
  USE mo_communication, ONLY: t_comm_pattern, delete_comm_pattern, &
    & exchange_data, exchange_data_mult
  USE mo_communication_factory, ONLY: setup_comm_pattern
  USE mo_fortran_tools, ONLY: set_acc_host_or_device
#ifndef NOMPI
  USE mpi, ONLY: MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE, MPI_COMM_NULL, MPI_UNDEFINED, &
    & MPI_Wait, MPI_Comm_split, MPI_Comm_rank
# ifndef NO_MPI_CHOICE_ARG
  ! Cray's MPI doesn't export MPI_Waitall.
  USE mpi, ONLY: MPI_Isend, MPI_Recv, MPI_Waitall
# endif
#endif

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
#ifdef NOMPI

    CALL finish(routine, "subset transfer does not make sense without MPI")
#else
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: gID_tmp
    INTEGER, ALLOCATABLE, DIMENSION(:) :: tmp_owners_cl, &
      & tmp_gIDs_sv, tmp_owners_sv, tmp_gIDs_sv_2, tmp_owners_sv_2, &
      & cl_ncor, cl_nbnd, cl_ntot, cl_rnk
    INTEGER :: ierr, n_cl, i_cl, iidx, iblk, rnk, i, req_cl(5), &
      & n_glb, n_core_sv, n_alloc_sv, nidx, n_core_cl, dum(3), dum2(3)
    TYPE(t_glb2loc_index_lookup) :: lookupTable
    CLASS(t_comm_pattern), POINTER :: cpat_sync_c
    LOGICAL :: done
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: tmp_owners_cl, tmp_gIDs_sv, tmp_owners_sv
!DIR$ ATTRIBUTES ALIGN : 64 :: tmp_gIDs_sv_2, tmp_owners_sv_2, cl_ncor
!DIR$ ATTRIBUTES ALIGN : 64 :: cl_ntot, cl_rnk, gID_tmp, cl_nbnd
#endif

    IF (this%is_init) RETURN
    IF (ltimer) THEN
      this%timer_init = new_timer("triv-t init")
      CALL timer_start(this%timer_init)
    END IF
    SELECT CASE(mode)
    CASE(solve_trans_scatter)
      this%solve_buddy = p_pe_work - MOD(p_pe_work, redfac)
      this%is_solver_pe = this%solve_buddy .EQ. p_pe_work
    CASE(solve_trans_compact)
      this%solve_buddy = p_pe_work / redfac
      this%is_solver_pe = p_pe_work .LT. (p_n_work + redfac - 1) / redfac
    CASE DEFAULT
      CALL finish(routine, "subset mode not recognized")
    END SELECT
    this%comm = MPI_COMM_NULL
    CALL MPI_Comm_split(p_comm_work, MERGE(1, MPI_UNDEFINED, this%is_solver_pe), &
      & 1, this%comm, ierr)
    rnk = -1
    IF (this%is_solver_pe) CALL MPI_Comm_rank(this%comm, rnk, ierr)
    this%is_leader_pe = rnk .EQ. 0
    this%nidx_l = nproma
    SELECT CASE(sync_type)
    CASE(solve_cell)
      this%nblk_a_l = patch_2D%alloc_cell_blocks
      this%nblk_l = patch_2d%cells%in_domain%end_block
      this%nidx_e_l = patch_2d%cells%in_domain%end_index
      this%glb_idx_loc => patch_2d%cells%decomp_info%glb_index
      cpat_sync_c => patch_2d%comm_pat_c
    CASE(solve_edge)
      this%nblk_a_l = SIZE(patch_2D%edges%decomp_info%owner_mask, 2)
      this%nblk_l = patch_2d%edges%in_domain%end_block
      this%nidx_e_l = patch_2d%edges%in_domain%end_index
      this%glb_idx_loc => patch_2d%edges%decomp_info%glb_index
      cpat_sync_c => patch_2d%comm_pat_e
    CASE(solve_vert)
      this%nblk_a_l = SIZE(patch_2D%verts%decomp_info%owner_mask, 2)
      this%nblk_l = patch_2d%verts%in_domain%end_block
      this%nidx_e_l = patch_2d%verts%in_domain%end_index
      this%glb_idx_loc => patch_2d%verts%decomp_info%glb_index
      cpat_sync_c => patch_2d%comm_pat_v
    CASE DEFAULT
      CALL finish(routine, "syncing scheme not recognized")
    END SELECT
    this%ngid_a_l = this%nidx_l * this%nblk_a_l
    n_alloc_sv = 0
    n_core_sv = 0
    CALL MPI_Isend(this%ngid_a_l, 1, p_int, this%solve_buddy, 1, &
      & p_comm_work, req_cl(1), ierr)
    IF (this%is_solver_pe) THEN
      n_cl = MIN(redfac, p_n_work - rnk * redfac)
      ALLOCATE(cl_ncor(n_cl), cl_nbnd(n_cl), cl_ntot(n_cl), cl_rnk(n_cl))
      cl_rnk(:) = (/(rnk * redfac + i - 1, i = 1, n_cl)/)
      nidx = 0
      DO i_cl = 1, n_cl
        CALL MPI_Recv(cl_ntot(i_cl), 1, p_int, cl_rnk(i_cl), 1, &
          & p_comm_work, MPI_STATUS_IGNORE, ierr)
        n_alloc_sv = n_alloc_sv + cl_ntot(i_cl)
        nidx = MAX(cl_ntot(i_cl), nidx)
      END DO
      ALLOCATE(tmp_gIDs_sv(nidx), tmp_owners_sv(nidx))
    END IF
    ALLOCATE(gID_tmp(this%nidx_l, this%nblk_a_l))
    this%nidx = nproma
    DO iblk = 1, this%nblk_l
      nidx = MERGE(this%nidx_l, this%nidx_e_l, this%nblk_l .NE. iblk)
      DO iidx = 1, nidx
        gID_tmp(iidx, iblk) = this%globalID_loc(iidx, iblk)
      END DO
      IF (nidx .NE. this%nidx_l) &
        & gID_tmp(nidx+1:this%nidx_l, iblk) = -1
    END DO
    DO iblk = this%nblk_l + 1, this%nblk_a_l
      gID_tmp(1:this%nidx_l, iblk) = -1
    ENDDO
    n_core_cl = (this%nblk_l - 1) * this%nidx_l + this%nidx_e_l
    n_glb = p_sum(n_core_cl, comm=p_comm_work)
    IF (my_process_is_mpi_parallel()) &
      & CALL exchange_data(p_pat=cpat_sync_c, lacc=i_am_accel_node, recv=gID_tmp)
    this%ngid_a_l = COUNT(gID_tmp(:,:) .NE. -1)
    NULLIFY(this%glb_idx_loc)
    ALLOCATE(this%glb_idx_loc(this%ngid_a_l))
    done = .FALSE.
    DO iblk = 1, this%nblk_a_l
      DO iidx = 1, this%nidx_l
        i = iidx + (iblk - 1) * this%nidx_l
        IF (i .GT.  this%ngid_a_l) THEN
          done = .TRUE.
          EXIT
        END IF
        this%glb_idx_loc(i) = gID_tmp(iidx, iblk)
      END DO
      IF (done) EXIT
    END DO
    DO iblk = 1, this%nblk_l
      nidx = MERGE(this%nidx_l, this%nidx_e_l, this%nblk_l .NE. iblk)
      gID_tmp(1:nidx, iblk) = p_pe_work
      IF (nidx .NE. this%nidx_l) gID_tmp(nidx+1:this%nidx_l, iblk) = -1
    END DO
    DO iblk = this%nblk_l + 1, this%nblk_a_l
      gID_tmp(1:this%nidx_l, iblk) = -1
    ENDDO
    IF (my_process_is_mpi_parallel()) &
      & CALL exchange_data(p_pat=cpat_sync_c, lacc=i_am_accel_node, recv=gID_tmp)
    ALLOCATE(tmp_owners_cl(this%ngid_a_l))
    done = .FALSE.
    DO iblk = 1, this%nblk_a_l
      DO iidx = 1, this%nidx_l
        i = iidx + (iblk - 1) * this%nidx_l
        IF (i .GT. this%ngid_a_l) THEN
          done = .TRUE.
          EXIT
        END IF
        tmp_owners_cl(i) = gID_tmp(iidx, iblk)
      END DO
      IF (done) EXIT
    END DO
    this%ngid_a_l = COUNT(tmp_owners_cl(1:this%ngid_a_l) .NE. -1)
    dum(1) = this%ngid_a_l
    dum(2) = COUNT(tmp_owners_cl(1:this%ngid_a_l) .EQ. p_pe_work)
    dum(3) = COUNT(tmp_owners_cl(1:this%ngid_a_l) .NE. p_pe_work &
           & .AND. tmp_owners_cl(1:this%ngid_a_l) .NE. -1)
    dum2(:) = dum(:)
    DEALLOCATE(gID_tmp)
    CALL MPI_Wait(req_cl(1), MPI_STATUS_IGNORE, ierr)
    CALL MPI_Isend(dum2, 3, p_int, this%solve_buddy, 11, p_comm_work, req_cl(1), ierr)
    CALL MPI_Isend(this%glb_idx_loc, this%ngid_a_l, p_int, this%solve_buddy, 4, p_comm_work, req_cl(2), ierr)
    CALL MPI_Isend(tmp_owners_cl, this%ngid_a_l, p_int, this%solve_buddy, 5, p_comm_work, req_cl(3), ierr)
    CALL MPI_Isend(this%glb_idx_loc, this%ngid_a_l, p_int, this%solve_buddy, 6, p_comm_work, req_cl(4), ierr)
    CALL MPI_Isend(tmp_owners_cl, this%ngid_a_l, p_int, this%solve_buddy, 7, p_comm_work, req_cl(5), ierr)
    IF (this%is_solver_pe) THEN
      DO i_cl = 1, n_cl
        CALL MPI_Recv(dum, 3, p_int, cl_rnk(i_cl), 11, &
          & p_comm_work, MPI_STATUS_IGNORE, ierr)
        cl_ntot(i_cl) = dum(1)
        cl_ncor(i_cl) = dum(2)
        cl_nbnd(i_cl) = dum(3)
      END DO
      i = SUM(cl_ncor(:)) + SUM(cl_nbnd(:))
      ALLOCATE(tmp_gIDs_sv_2(i), tmp_owners_sv_2(i))
      DO i_cl = 1, n_cl
        CALL MPI_Recv(tmp_gIDs_sv, cl_ntot(i_cl), p_int, &
          & cl_rnk(i_cl), 4, p_comm_work, MPI_STATUS_IGNORE, ierr)
        CALL MPI_Recv(tmp_owners_sv, cl_ntot(i_cl), p_int, &
          & cl_rnk(i_cl), 5, p_comm_work, MPI_STATUS_IGNORE, ierr)
        DO i = 1, cl_ncor(i_cl)
          tmp_gIDs_sv_2(n_core_sv + i) = tmp_gIDs_sv(i)
          tmp_owners_sv_2(n_core_sv + i) = cl_rnk(i_cl)
        END DO
        n_core_sv = n_core_sv + cl_ncor(i_cl)
      END DO
      n_alloc_sv = n_core_sv
      DO i_cl = 1, n_cl
        CALL MPI_Recv(tmp_gIDs_sv, cl_ntot(i_cl), p_int, &
          & cl_rnk(i_cl), 6, p_comm_work, MPI_STATUS_IGNORE, ierr)
        CALL MPI_Recv(tmp_owners_sv, cl_ntot(i_cl), p_int, &
          & cl_rnk(i_cl), 7, p_comm_work, MPI_STATUS_IGNORE, ierr)
        DO i = cl_ncor(i_cl) + 1, cl_ncor(i_cl) + cl_nbnd(i_cl)
          IF (tmp_owners_sv(i) .NE. -1) THEN
            IF (.NOT.ANY(tmp_gIDs_sv_2(1:n_alloc_sv) .EQ. tmp_gIDs_sv(i))) THEN
              n_alloc_sv = n_alloc_sv + 1
              tmp_gIDs_sv_2(n_alloc_sv) = tmp_gIDs_sv(i)
              tmp_owners_sv_2(n_alloc_sv) = tmp_owners_sv(i)
            END IF
          END IF
        END DO
      END DO
      DEALLOCATE(tmp_owners_sv, tmp_gIDs_sv)
      ALLOCATE(this%glb_idx_cal(n_alloc_sv), tmp_owners_sv(n_alloc_sv))
      DO i = 1, n_alloc_sv
        this%glb_idx_cal(i) = tmp_gIDs_sv_2(i)
      END DO
      DEALLOCATE(tmp_gIDs_sv_2, cl_ncor, cl_nbnd, cl_ntot)
    ELSE
      ALLOCATE(this%glb_idx_cal(0), tmp_owners_sv(0), tmp_owners_sv_2(0))
    END IF
    IF (this%solve_buddy .NE. p_pe_work) &
      & CALL MPI_Waitall(5, req_cl, MPI_STATUSES_IGNORE, ierr)
    this%ngid_a = n_alloc_sv
    CALL init_glb2loc_index_lookup(lookupTable, n_glb)
    CALL set_inner_glb_index(lookupTable, this%glb_idx_loc(1:n_core_cl), &
      & [(i, i = 1, n_core_cl)])
    CALL setup_comm_pattern(n_alloc_sv, tmp_owners_sv_2, this%glb_idx_cal, lookupTable, &
      & n_core_cl, tmp_owners_cl, this%glb_idx_loc, this%cpat_in, &
      & inplace=.FALSE., comm=p_comm_work)
    CALL setup_comm_pattern(n_core_sv, tmp_owners_sv_2(1:n_core_sv), &
      & this%glb_idx_cal(1:n_core_sv), lookupTable, n_core_cl, tmp_owners_cl(1:n_core_cl), &
      & this%glb_idx_loc(1:n_core_cl), this%cpat_in2, inplace=.FALSE., comm=p_comm_work)
    CALL deallocate_glb2loc_index_lookup(lookupTable)
    SELECT CASE(mode)
    CASE(solve_trans_scatter)
      tmp_owners_cl(:) = MERGE(-1, tmp_owners_cl(:) - MOD(tmp_owners_cl(:), redfac), &
        & tmp_owners_cl(:) .EQ. -1)
      IF (this%is_solver_pe) &
        & tmp_owners_sv(1:n_core_sv) = MERGE(-1, &
          & tmp_owners_sv_2(1:n_core_sv) - MOD(tmp_owners_sv_2(1:n_core_sv), redfac), &
          & tmp_owners_sv_2(1:n_core_sv) .EQ. -1)
    CASE(solve_trans_compact)
      tmp_owners_cl(:) = MERGE(-1, tmp_owners_cl(:) / redfac, tmp_owners_cl(:) .EQ. -1)
      IF (this%is_solver_pe) &
        & tmp_owners_sv(1:n_core_sv) = MERGE(-1, tmp_owners_sv_2(1:n_core_sv) / redfac, &
          & tmp_owners_sv_2(1:n_core_sv) .EQ. -1)
    END SELECT
    CALL init_glb2loc_index_lookup(lookupTable, n_glb)
    CALL set_inner_glb_index(lookupTable, this%glb_idx_cal(1:n_core_sv), &
      & [(i, i = 1, n_core_sv)])
    CALL setup_comm_pattern(this%ngid_a_l, tmp_owners_cl, this%glb_idx_loc, &
        & lookupTable, n_core_sv, tmp_owners_sv(1:n_core_sv), &
        & this%glb_idx_cal(1:n_core_sv), this%cpat_out, inplace=.FALSE., comm=p_comm_work)
    IF (this%is_solver_pe) THEN
      tmp_owners_sv(1:n_core_sv) = rnk
      tmp_owners_sv(n_core_sv + 1:) = -1
      tmp_owners_sv_2(1:n_core_sv) = -1
      tmp_owners_sv_2(n_core_sv + 1:) = tmp_owners_sv_2(n_core_sv + 1:) / redfac
      CALL setup_comm_pattern(n_alloc_sv, tmp_owners_sv_2, this%glb_idx_cal, lookupTable, &
          & n_core_sv, tmp_owners_sv(1:n_core_sv), &
          & this%glb_idx_cal(1:n_core_sv), this%cpat_sync, inplace=.TRUE., comm=this%comm)
    END IF
    CALL deallocate_glb2loc_index_lookup(lookupTable)
    DEALLOCATE(tmp_owners_sv, tmp_owners_sv_2)
    this%nblk_a = (n_alloc_sv + nproma - 1) / nproma
    this%nblk = (n_core_sv + nproma - 1) / nproma
    this%nidx_e = n_core_sv - (this%nblk - 1) * nproma
    IF (ltimer) THEN
      this%timer_in(1) = new_timer("sub-t in(solve)")
      this%timer_in(2) = new_timer("sub-t in(lhs)")
      this%timer_in(3) = new_timer("sub-t in(init solve)")
      this%timer_in(4) = new_timer("sub-t in(init lhs)")
      this%timer_out = new_timer("sub-t out")
      CALL timer_stop(this%timer_init)
      this%timer_glob_sum = new_timer("sub-t glb_sum")
      this%timer_sync = new_timer("sub-t sync")
    END IF
    this%is_init = .true.
#endif !! ifndef NOMPI
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

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

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

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

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

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

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

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

    IF (ltimer) CALL timer_start(this%timer_in(tt))
    CALL exchange_data(p_pat=this%cpat_in, lacc=lzacc, recv=data_out, send=data_in)
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

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: to, ti
#endif

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

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
    CALL exchange_data_mult(p_pat=this%cpat_in, lacc=lzacc, nfields=2, ndim2tot=2, recv4d=to, send4d=ti)
    IF (this%is_solver_pe) THEN
#ifdef _OPENMP
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
#endif
      DO i = 1, 2
        IF (i .EQ. 1) THEN
          do1(:,:) = to(:,1,:,i)
        ELSE
          do2(:,:) = to(:,1,:,i)
        END IF
      END DO
#ifdef _OPENMP
!ICON_OMP END PARALLEL DO
#endif
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

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: to
#endif

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

    IF (ltimer) CALL timer_start(this%timer_in(tt))
    n3 = SIZE(data_in, 3)
    ALLOCATE(to(this%nidx, 1, this%nblk, n3))
    IF (this%is_solver_pe) to(:, 1, this%nblk, :) = 0._wp
    CALL exchange_data_mult(p_pat=this%cpat_in2, lacc=lzacc, nfields=n3, ndim2tot=n3, recv4d=to, send4d= &
      & RESHAPE(data_in, (/SIZE(data_in, 1), 1, SIZE(data_in, 2), n3/)))
    IF (this%is_solver_pe) THEN
#ifdef _OPENMP
!$OMP PARALLEL DO SCHEDULE(STATIC)
#endif
      DO i = 1, n3
        data_out(:,:,i) = to(:,1,:,i)
      END DO
#ifdef _OPENMP
!$OMP END PARALLEL DO
#endif
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

#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: glb_in, glb_out
#endif

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

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
      CALL exchange_data(p_pat=this%cpat_in2, lacc=lzacc, recv=glb_out, send=glb_in)
      IF (this%is_solver_pe) THEN
#ifdef _OPENMP
!$OMP PARALLEL DO SCHEDULE(GUIDED) PRIVATE(iidx, jblk, jidx, gid, found, notfound)
#endif
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
#ifdef _OPENMP
!$OMP END PARALLEL DO
#endif
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

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

    IF (ltimer) CALL timer_start(this%timer_out)
    CALL exchange_data(p_pat=this%cpat_out, lacc=lzacc, recv=data_out, send=data_in)
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

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

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

#ifdef _OPENACC
    if (lzacc) CALL finish(routine, "not ported to GPU yet")
#endif

    IF (ltimer) CALL timer_start(this%timer_out)
    IF (this%is_leader_pe) data_out(:) = data_in(:)
    CALL p_bcast(data_out, 0)
    IF (ltimer) CALL timer_stop(this%timer_out)
  END SUBROUTINE subset_transfer_bcst_1d_i

  SUBROUTINE subset_transfer_sync_2d_wp(this, data_inout)
    CLASS(t_subset_transfer), INTENT(INOUT) :: this
    REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout

    IF (ltimer) CALL timer_start(this%timer_sync)
    IF (my_process_is_mpi_parallel()) THEN
#ifdef _OPENACC
      CALL finish("subset_transfer_sync_2d_wp is not tested with OPENACC")
      ! TODO: lacc=.flase. or .true.?
#endif
      CALL exchange_data(p_pat=this%cpat_sync, lacc=i_am_accel_node, recv=data_inout)
    END IF
    IF (ltimer) CALL timer_stop(this%timer_sync)
  END SUBROUTINE subset_transfer_sync_2d_wp

  SUBROUTINE subset_transfer_sync_2d_sp(this, data_inout)
    CLASS(t_subset_transfer), INTENT(INOUT) :: this
    REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout

    IF (ltimer) CALL timer_start(this%timer_sync)
    IF (my_process_is_mpi_parallel()) THEN
#ifdef _OPENACC
      CALL finish("subset_transfer_sync_2d_sp is not tested with OPENACC")
      ! TODO: lacc=.flase. or .true.?
#endif
      CALL exchange_data(p_pat=this%cpat_sync, lacc=i_am_accel_node, recv=data_inout)
    END IF
    IF (ltimer) CALL timer_stop(this%timer_sync)
  END SUBROUTINE subset_transfer_sync_2d_sp

END MODULE mo_ocean_solve_subset_transfer
