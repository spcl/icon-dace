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

! provides abstract communication / transfer infrastructure object to
! be used by solvers

#if (defined(_OPENMP) && defined(OCE_SOLVE_OMP))
#include "omp_definitions.inc"
#define PURE_OR_OMP
#else
#define PURE_OR_OMP PURE
#endif

MODULE mo_ocean_solve_transfer

  USE mo_kind, ONLY: sp, dp, i8
  USE mo_fortran_tools, ONLY: t_ptr_2d, t_ptr_2d_sp
  USE mo_parallel_config, ONLY: l_fast_sum
  USE mo_timer, ONLY: timer_start, timer_stop
  USE mo_mpi, ONLY: p_sum, p_max
  USE mo_exception, ONLY: finish
  USE mo_run_config, ONLY: ltimer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_transfer

  TYPE, ABSTRACT :: t_transfer
    PRIVATE
    INTEGER, PUBLIC :: comm, nblk_a, nblk, nidx, nidx_e, nidx_l, ngid_a_l
    INTEGER, PUBLIC :: timer_sync, timer_in(4), timer_out, timer_glob_sum, timer_init
    LOGICAL, PUBLIC :: is_solver_pe, is_leader_pe
    LOGICAL, PUBLIC :: is_init = .false.
    INTEGER, POINTER, DIMENSION(:), PUBLIC :: glb_idx_loc => NULL(), &
      & glb_idx_cal => NULL()
  CONTAINS
! constructor method must be defined in extended type ... not declared deferred,
! since different implementations may use different constructor args
! transfer array data from worker- to solver-PEs, a gather operation, if n_solverPE < n_workerPE
! (assumes arrays on solverside are already allocated)
    PROCEDURE(a_trans_into_2d_wp), DEFERRED :: into_2d_wp
    PROCEDURE(a_trans_into_2d_wp_2), DEFERRED :: into_2d_wp_2
    PROCEDURE(a_trans_into_3d_wp), DEFERRED :: into_3d_wp
    PROCEDURE(a_trans_into_idx), DEFERRED :: into_idx
    GENERIC, PUBLIC :: into => into_2d_wp, into_2d_wp_2, into_3d_wp, into_idx
! transfer array data from worker- to solver-PEs (first transfer -> allocates arrays on solver-PE side)
    PROCEDURE(a_trans_into_once_2d_wp), DEFERRED :: into_once_2d_wp
    PROCEDURE(a_trans_into_once_3d_wp), DEFERRED :: into_once_3d_wp
    PROCEDURE(a_trans_into_once_idx), DEFERRED :: into_once_idx
    GENERIC, PUBLIC :: into_once => into_once_2d_wp, into_once_3d_wp, into_once_idx
! transfer array data back from solve- to worker-PEs (a scatter operation)
    PROCEDURE(a_trans_out_2d_wp), DEFERRED :: out_2d_wp
    GENERIC, PUBLIC :: sctr => out_2d_wp
! broadcast array data from solve- to worker-PEs (e.g. for residual-norms)
    PROCEDURE(a_trans_bcst_1d_wp), DEFERRED :: bcst_1d_wp
    PROCEDURE(a_trans_bcst_1d_i), DEFERRED :: bcst_1d_i
    GENERIC, PUBLIC :: bcst => bcst_1d_wp, bcst_1d_i
! boundary sync between solver-PEs
    PROCEDURE(a_trans_sync_2d_wp), DEFERRED :: sync_2d_wp
    PROCEDURE(a_trans_sync_2d_sp), DEFERRED :: sync_2d_sp
    GENERIC, PUBLIC :: sync => sync_2d_wp, sync_2d_sp
! global sums (for sp + wp real-kinds, and for up to 3 sums at once)
    PROCEDURE, PRIVATE :: gs_2d_1dp => ocean_solve_transfer_global_sum_2d_dp_1
    PROCEDURE, PRIVATE :: gs_2d_2dp => ocean_solve_transfer_global_sum_2d_dp_2
    PROCEDURE, PRIVATE :: gs_2d_3dp => ocean_solve_transfer_global_sum_2d_dp_3
    PROCEDURE, PRIVATE :: gs_2d_1sp => ocean_solve_transfer_global_sum_2d_sp_1
    PROCEDURE, PRIVATE :: gs_2d_2sp => ocean_solve_transfer_global_sum_2d_sp_2
    PROCEDURE, PRIVATE :: gs_2d_3sp => ocean_solve_transfer_global_sum_2d_sp_3
    GENERIC, PUBLIC :: global_sum => gs_2d_1dp, gs_2d_1sp, &
      & gs_2d_2dp, gs_2d_2sp, gs_2d_3dp, gs_2d_3sp
! global sums internal routines (for sp + wp real-kinds, and for N sums at once)
    PROCEDURE, PRIVATE :: gs_2d_ndp => ocean_solve_transfer_global_sum_2d_dp
    PROCEDURE, PRIVATE :: gs_2d_nsp => ocean_solve_transfer_global_sum_2d_sp
    GENERIC, PRIVATE :: global_sum_internal => gs_2d_ndp, gs_2d_nsp
! global indices
    PROCEDURE :: globalID_loc => ocean_solve_transfer_globalID_loc
    PROCEDURE :: globalID_cal => ocean_solve_transfer_globalID_cal
    PROCEDURE(a_trans_destruct), DEFERRED :: destruct
  END TYPE t_transfer

  ABSTRACT INTERFACE
    SUBROUTINE a_trans_into_once_2d_wp(this, data_in, data_out, tt, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: data_in
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), ALLOCATABLE :: data_out
      INTEGER, INTENT(IN) :: tt
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_into_once_2d_wp
    SUBROUTINE a_trans_into_2d_wp_2(this, di1, do1, di2, do2, tt, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: di1, di2
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: do1, do2
      INTEGER, INTENT(IN) :: tt
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_into_2d_wp_2
    SUBROUTINE a_trans_into_once_3d_wp(this, data_in, data_out, tt, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:,:), ALLOCATABLE :: data_out
      INTEGER, INTENT(IN) :: tt
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_into_once_3d_wp
    SUBROUTINE a_trans_into_once_idx(this, data_in_idx, data_in_blk, &
      & data_out_idx, data_out_blk, tt, lacc)
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in_idx, data_in_blk
      INTEGER, INTENT(INOUT), DIMENSION(:,:,:), ALLOCATABLE :: &
        & data_out_idx, data_out_blk
      INTEGER, INTENT(IN) :: tt
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_into_once_idx
    SUBROUTINE a_trans_into_2d_wp(this, data_in, data_out, tt, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      REAL(KIND=wp), INTENT(IN), DIMENSION(:,:) :: data_in
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_out
      INTEGER, INTENT(IN) :: tt
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_into_2d_wp
    SUBROUTINE a_trans_into_3d_wp(this, data_in, data_out, tt, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      REAL(KIND=wp), INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:,:), CONTIGUOUS :: data_out
      INTEGER, INTENT(IN) :: tt
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_into_3d_wp
    SUBROUTINE a_trans_into_idx(this, data_in_idx, data_in_blk, &
      & data_out_idx, data_out_blk, tt, lacc)
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      INTEGER, INTENT(IN), DIMENSION(:,:,:), CONTIGUOUS :: data_in_blk, data_in_idx
      INTEGER, INTENT(INOUT), DIMENSION(:,:,:), CONTIGUOUS :: data_out_blk, data_out_idx
      INTEGER, INTENT(IN) :: tt
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_into_idx
    SUBROUTINE a_trans_out_2d_wp(this, data_in, data_out, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      REAL(KIND=wp), INTENT(IN), DIMENSION(:,:), CONTIGUOUS :: data_in
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_out
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_out_2d_wp
    SUBROUTINE a_trans_bcst_1d_wp(this, data_in, data_out, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      REAL(KIND=wp), INTENT(IN), DIMENSION(:), CONTIGUOUS :: data_in
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: data_out
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_bcst_1d_wp
    SUBROUTINE a_trans_bcst_1d_i(this, data_in, data_out, lacc)
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(IN) :: this
      INTEGER, INTENT(IN), DIMENSION(:), CONTIGUOUS :: data_in
      INTEGER, INTENT(INOUT), DIMENSION(:), CONTIGUOUS :: data_out
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_bcst_1d_i
    SUBROUTINE a_trans_sync_2d_wp(this, data_inout, lacc)
      USE mo_kind, ONLY: wp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(INOUT) :: this
      REAL(KIND=wp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_sync_2d_wp
    SUBROUTINE a_trans_sync_2d_sp(this, data_inout, lacc)
      USE mo_kind, ONLY: sp
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(INOUT) :: this
      REAL(KIND=sp), INTENT(INOUT), DIMENSION(:,:), CONTIGUOUS :: data_inout
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_sync_2d_sp
    SUBROUTINE a_trans_destruct(this, lacc)
      IMPORT t_transfer
      CLASS(t_transfer), INTENT(INOUT) :: this
      LOGICAL, INTENT(IN), OPTIONAL :: lacc
    END SUBROUTINE a_trans_destruct
  END INTERFACE

  INTERFACE simple_sum_local
    MODULE PROCEDURE simple_sum_loc_dp_2d
    MODULE PROCEDURE simple_sum_loc_sp_2d
  END INTERFACE simple_sum_local

  INTERFACE abs_max_loc
    MODULE PROCEDURE abs_max_loc_dp_2d
    MODULE PROCEDURE abs_max_loc_sp_2d
  END INTERFACE abs_max_loc

  INTERFACE order_insensit_ieee64_sum_frst
    MODULE PROCEDURE order_insensit_ieee64_sum_frst_dp_2d
    MODULE PROCEDURE order_insensit_ieee64_sum_frst_sp_2d
  END INTERFACE order_insensit_ieee64_sum_frst

  INTERFACE order_insensit_ieee64_sum_scnd
    MODULE PROCEDURE order_insensit_ieee64_sum_scnd_2d
  END INTERFACE order_insensit_ieee64_sum_scnd

  CHARACTER(LEN=*), PARAMETER :: module_name = "mo_ocean_solve_transfer"

CONTAINS

! returns global index for (iidx,iblk) in a worker-side array
  PURE ELEMENTAL FUNCTION ocean_solve_transfer_globalID_loc(this, &
      & idx, blk) RESULT(gidx)
    CLASS(t_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: idx, blk
    INTEGER :: gidx, lidx

    lidx = (blk-1)*this%nidx_l + idx
    gidx = this%glb_idx_loc(lidx)
  END FUNCTION ocean_solve_transfer_globalID_loc

! returns global index for (iidx,iblk) in a solver-side array
  PURE ELEMENTAL FUNCTION ocean_solve_transfer_globalID_cal(this, &
      & idx, blk) RESULT(gidx)
    CLASS(t_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: idx, blk
    INTEGER :: gidx, lidx

    lidx = (blk-1)*this%nidx + idx
    gidx = this%glb_idx_cal(lidx)
  END FUNCTION ocean_solve_transfer_globalID_cal

! explicite interface for 1 global sum of double precision
  SUBROUTINE ocean_solve_transfer_global_sum_2d_dp_1(this, &
    & data_in1, summa1)
    CLASS(t_transfer), INTENT(IN) :: this
    REAL(KIND=dp), INTENT(IN), TARGET :: data_in1(:,:)
    REAL(KIND=dp), INTENT(OUT) :: summa1
    TYPE(t_ptr_2d) :: data_in_ptr(1)
    REAL(KIND=dp) :: sums(1)
#ifdef __INTEL_COMPILER
    !DIR$ ATTRIBUTES ALIGN : 64 :: sums
#endif

    data_in_ptr(1)%p => data_in1
    CALL this%global_sum_internal(1, data_in_ptr, sums)
    summa1 = sums(1)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_dp_1

! explicite interface for 2 global sums of double precision
  SUBROUTINE ocean_solve_transfer_global_sum_2d_dp_2(this, &
    & data_in1, summa1, data_in2, summa2)
    CLASS(t_transfer), INTENT(IN) :: this
    REAL(KIND=dp), INTENT(IN), TARGET :: &
      & data_in1(:,:), data_in2(:,:)
    REAL(KIND=dp), INTENT(OUT) :: summa1, summa2
    TYPE(t_ptr_2d) :: data_in_ptr(2)
    REAL(KIND=dp) :: sums(2)
#ifdef __INTEL_COMPILER
    !DIR$ ATTRIBUTES ALIGN : 64 :: sums
#endif

    data_in_ptr(1)%p => data_in1
    data_in_ptr(2)%p => data_in2
    CALL this%global_sum_internal(2, data_in_ptr, sums)
    summa1 = sums(1)
    summa2 = sums(2)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_dp_2

! explicite interface for 3 global sums of double precision
  SUBROUTINE ocean_solve_transfer_global_sum_2d_dp_3(this, &
    & data_in1, summa1, data_in2, summa2, data_in3, summa3)
    CLASS(t_transfer), INTENT(IN) :: this
    REAL(KIND=dp), INTENT(IN), TARGET :: &
      & data_in1(:,:), data_in2(:,:), data_in3(:,:)
    REAL(KIND=dp), INTENT(OUT) :: summa1, summa2, summa3
    TYPE(t_ptr_2d) :: data_in_ptr(3)
    REAL(KIND=dp) :: sums(3)
#ifdef __INTEL_COMPILER
    !DIR$ ATTRIBUTES ALIGN : 64 :: sums
#endif

    data_in_ptr(1)%p => data_in1
    data_in_ptr(2)%p => data_in2
    data_in_ptr(3)%p => data_in3
    CALL this%global_sum_internal(3, data_in_ptr, sums)
    summa1 = sums(1)
    summa2 = sums(2)
    summa3 = sums(3)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_dp_3

! explicite interface for 1 global sum of single precision
  SUBROUTINE ocean_solve_transfer_global_sum_2d_sp_1(this, &
    & data_in1, summa1)
    CLASS(t_transfer), INTENT(IN) :: this
    REAL(KIND=sp), INTENT(IN), TARGET :: data_in1(:,:)
    REAL(KIND=sp), INTENT(OUT) :: summa1
    TYPE(t_ptr_2d_sp) :: data_in_ptr(1)
    REAL(KIND=sp) :: sums(1)
#ifdef __INTEL_COMPILER
    !DIR$ ATTRIBUTES ALIGN : 64 :: sums
#endif

    data_in_ptr(1)%p => data_in1
    CALL this%global_sum_internal(1, data_in_ptr, sums)
    summa1 = sums(1)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_sp_1

! explicite interface for 2 global sums of single precision
  SUBROUTINE ocean_solve_transfer_global_sum_2d_sp_2(this, &
    & data_in1, summa1, data_in2, summa2)
    CLASS(t_transfer), INTENT(IN) :: this
    REAL(KIND=sp), INTENT(IN), TARGET :: &
      & data_in1(:,:), data_in2(:,:)
    REAL(KIND=sp), INTENT(OUT) :: summa1, summa2
    TYPE(t_ptr_2d_sp) :: data_in_ptr(2)
    REAL(KIND=sp) :: sums(2)
#ifdef __INTEL_COMPILER
    !DIR$ ATTRIBUTES ALIGN : 64 :: sums
#endif

    data_in_ptr(1)%p => data_in1
    data_in_ptr(2)%p => data_in2
    CALL this%global_sum_internal(2, data_in_ptr, sums)
    summa1 = sums(1)
    summa2 = sums(2)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_sp_2

! explicite interface for 3 global sums of single precision
  SUBROUTINE ocean_solve_transfer_global_sum_2d_sp_3(this, &
    & data_in1, summa1, data_in2, summa2, data_in3, summa3)
    CLASS(t_transfer), INTENT(IN) :: this
    REAL(KIND=sp), INTENT(IN), TARGET :: &
      & data_in1(:,:), data_in2(:,:), data_in3(:,:)
    REAL(KIND=sp), INTENT(OUT) :: summa1, summa2, summa3
    TYPE(t_ptr_2d_sp) :: data_in_ptr(3)
    REAL(KIND=sp) :: sums(3)
#ifdef __INTEL_COMPILER
    !DIR$ ATTRIBUTES ALIGN : 64 :: sums
#endif

    data_in_ptr(1)%p => data_in1
    data_in_ptr(2)%p => data_in2
    data_in_ptr(3)%p => data_in3
    CALL this%global_sum_internal(3, data_in_ptr, sums)
    summa1 = sums(1)
    summa2 = sums(2)
    summa3 = sums(3)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_sp_3

! explicite internal interface for N global sums of double precision
! uses order insensitive or 'fast' implementation depending on l_fast_sum
  SUBROUTINE ocean_solve_transfer_global_sum_2d_dp(this, n, xp, gbl_sum)
    CLASS(t_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: n
    TYPE(t_ptr_2d), INTENT(IN) :: xp(n)
    REAL(dp), INTENT(OUT) :: gbl_sum(n)
    REAL(dp) :: loc_sum(n), abs_max_l(n), abs_max(n)
    INTEGER(KIND=i8) :: isum_loc(2*n), isum(2*n), tisum(2)
    INTEGER :: i
    CHARACTER(LEN=*), PARAMETER :: routine = module_name // &
      & '::global_sum_2d_dp'
#ifdef __INTEL_COMPILER
    !DIR$ ATTRIBUTES ALIGN : 64 :: loc_sum, abs_max_l, abs_max
#endif
    gbl_sum(:) = 0._dp
    loc_sum(:) = 0._dp
    IF (ltimer) CALL timer_start(this%timer_glob_sum)
    IF (l_fast_sum) THEN
      DO i = 1, n
#ifdef _CRAYFTN
!DIR$ NOINLINE
#endif
        loc_sum(i) = simple_sum_local(xp(i)%p)
      END DO
      gbl_sum(:) = p_sum(loc_sum(:), comm=this%comm)
    ELSE
      DO i = 1, n
        abs_max_l(i) = abs_max_loc(xp(i)%p)
      END DO
      abs_max(:) = p_max(abs_max_l(:), comm=this%comm)
      DO i = 1, n
        tisum(:) = order_insensit_ieee64_sum_frst( &
          & xp(i)%p, abs_max(i))
        isum_loc((i-1)*2+1:i*2) = tisum(:)
      END DO
      isum(:) = p_sum(isum_loc(:), comm=this%comm)
      DO i = 1, n
        tisum(:) = isum((i-1)*2+1:i*2)
        gbl_sum(i) = order_insensit_ieee64_sum_scnd( &
          & tisum(:), abs_max(i))
      END DO
    END IF
    IF (ltimer) CALL timer_stop(this%timer_glob_sum)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_dp

! explicite internal interface for N global sums of single precision
! uses order insensitive or 'fast' implementation depending on l_fast_sum
  SUBROUTINE ocean_solve_transfer_global_sum_2d_sp(this, n, xp, gbl_sum)
    CLASS(t_transfer), INTENT(IN) :: this
    INTEGER, INTENT(IN) :: n
    TYPE(t_ptr_2d_sp), INTENT(IN) :: xp(n)
    REAL(sp), INTENT(OUT) :: gbl_sum(n)
    REAL(sp) :: loc_sum(n), abs_max_l(n), abs_max(n)
    INTEGER :: i
    INTEGER(KIND=i8) :: isum_loc(2*n), isum(2*n), tisum(2)
    CHARACTER(LEN=*), PARAMETER :: routine = module_name // &
      & '::global_sum_2d_sp'
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: loc_sum, abs_max_l, abs_max
#endif
    gbl_sum(:) = 0._sp
    loc_sum(:) = 0._sp
    IF (ltimer) CALL timer_start(this%timer_glob_sum)
    IF(l_fast_sum) THEN
      DO i = 1, n
#ifdef _CRAYFTN
!DIR$ NOINLINE
#endif
        loc_sum(i) = simple_sum_local(xp(i)%p)
      END DO
      gbl_sum(:) = p_sum(loc_sum(:), comm=this%comm)
    ELSE
      DO i = 1, n
        abs_max_l(i) = abs_max_loc(xp(i)%p)
      END DO
      abs_max(:) = p_max(abs_max_l(:), comm=this%comm)
      DO i = 1, n
        tisum(:) = order_insensit_ieee64_sum_frst( &
          & xp(i)%p, abs_max(i))
        isum_loc((i-1)*2+1:i*2) = tisum(:)
      END DO
      isum(:) = p_sum(isum_loc(:), comm=this%comm)
      DO i = 1, n
        tisum(:) = isum((i-1)*2+1:i*2)
        gbl_sum(i) = REAL(order_insensit_ieee64_sum_scnd( &
          & tisum(:), REAL(abs_max(i), dp)), sp)
      END DO
    END IF
    IF (ltimer) CALL timer_stop(this%timer_glob_sum)
  END SUBROUTINE ocean_solve_transfer_global_sum_2d_sp

! performs local sum -- 'fast' implementation - dp variant
  PURE_OR_OMP FUNCTION simple_sum_loc_dp_2d(vals) RESULT(local_sum)
    REAL(dp), INTENT(IN), CONTIGUOUS :: vals(:,:)
    REAL(dp) :: local_sum, aux_sum(SIZE(vals, 2))
    INTEGER :: j
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: aux_sum
#endif

!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO j = 1, SIZE(vals, 2)
      aux_sum(j) = SUM(vals(:,j))
    END DO
    local_sum = SUM(aux_sum(:))
  END FUNCTION simple_sum_loc_dp_2d

! performs local sum -- 'fast' implementation - sp variant
  PURE_OR_OMP FUNCTION simple_sum_loc_sp_2d(vals) RESULT(local_sum)
    REAL(sp), INTENT(IN), CONTIGUOUS :: vals(:,:)
    REAL(sp) :: local_sum, aux_sum(SIZE(vals, 2))
    INTEGER :: j
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: aux_sum
#endif

!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO j = 1, SIZE(vals, 2)
      aux_sum(j) = SUM(vals(:,j))
    END DO
    local_sum = SUM(aux_sum(:))
  END FUNCTION simple_sum_loc_sp_2d

! finds local MAXVAL(ABS(x(:,:)) implementation - dp variant
  PURE_OR_OMP FUNCTION abs_max_loc_dp_2d(vals) RESULT(abs_max)
    REAL(dp), INTENT(IN), CONTIGUOUS :: vals(:,:)
    REAL(dp) :: abs_max, aux_max(SIZE(vals, 2))
    INTEGER :: j
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: aux_max
#endif

!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO j = 1, SIZE(vals, 2)
      aux_max(j) = MAXVAL(ABS(vals(:,j)))
    END DO
    abs_max = MAXVAL(aux_max(:))
  END FUNCTION abs_max_loc_dp_2d

! finds local MAXVAL(ABS(x(:,:)) implementation - dp variant
  PURE_OR_OMP FUNCTION abs_max_loc_sp_2d(vals) RESULT(abs_max)
    REAL(sp), INTENT(IN), CONTIGUOUS :: vals(:,:)
    REAL(sp) :: abs_max, aux_max(SIZE(vals, 2))
    INTEGER :: j
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: aux_max
#endif

!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO j = 1, SIZE(vals, 2)
      aux_max(j) = MAXVAL(ABS(vals(:,j)))
    END DO
    abs_max = MAXVAL(aux_max(:))  
END FUNCTION abs_max_loc_sp_2d

! first local part of order insensitive summation -- dp-variant
! convert to scaled integers and locally sum those
  PURE_OR_OMP FUNCTION order_insensit_ieee64_sum_frst_dp_2d(vals, abs_max) &
    & RESULT(isum)
    REAL(dp), INTENT(IN), CONTIGUOUS :: vals(:,:)
    REAL(dp), INTENT(IN) :: abs_max
    INTEGER(KIND=i8) :: isum(2), isum1(SIZE(vals, 2)), &
      & isum2(SIZE(vals, 2)), ival(SIZE(vals, 1))
    INTEGER :: j, iexp
    REAL(KIND=dp) :: fact, rval(SIZE(vals, 1))
    REAL(dp), PARAMETER :: two_30 = 1073741824._dp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: rval, ival, isum1, isum2
#endif

    iexp = EXPONENT(abs_max)
    IF (iexp < -980) THEN
      isum(:) = 0_i8
      RETURN
    END IF
    fact = SCALE(1._dp,30-iexp)
    isum1(:) = 0_i8
    isum2(:) = 0_i8
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ival, rval)
    DO j = 1, SIZE(vals, 2)
      rval(:) = vals(:,j) * fact
      ival(:) = INT(rval(:), i8)
      isum1(j) = isum1(j) + SUM(ival(:))
      isum2(j) = isum2(j) + SUM(INT((rval(:) - REAL(ival(:),dp))*two_30,i8))
    END DO
    isum(1) = SUM(isum1(:))
    isum(2) = SUM(isum2(:))
  END FUNCTION order_insensit_ieee64_sum_frst_dp_2d

! first local part of order insensitive summation -- sp-variant
! convert to scaled integers and locally sum those
  FUNCTION order_insensit_ieee64_sum_frst_sp_2d(vals, abs_max) &
    & RESULT(isum)
    REAL(sp), INTENT(IN), CONTIGUOUS :: vals(:,:)
    REAL(sp), INTENT(IN) :: abs_max
    INTEGER(KIND=i8) :: isum(2), isum1(SIZE(vals, 2)), &
      & isum2(SIZE(vals, 2)), ival(SIZE(vals, 1)) 
    INTEGER :: j, iexp
    REAL(KIND=dp) :: fact, rval(SIZE(vals, 1))
    REAL(dp), PARAMETER :: two_30 = 1073741824._dp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: rval, ival, isum1, isum2
#endif

    iexp = EXPONENT(abs_max)
    IF (iexp < -980) THEN
      isum(:) = 0_i8
      RETURN
    END IF
    fact = SCALE(1._dp,30-iexp)
    isum1(:) = 0_i8
    isum2(:) = 0_i8
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ival, rval)
    DO j = 1, SIZE(vals, 2)
      rval(:) = REAL(vals(:,j), dp) * fact
      ival(:) = INT(rval(:), i8)
      isum1(j) = isum1(j) + SUM(ival(:))
      isum2(j) = isum2(j) + SUM(INT((rval(:) - REAL(ival(:),dp))*two_30,i8))
    END DO
    isum(1) = SUM(isum1(:))
    isum(2) = SUM(isum2(:))
  END FUNCTION order_insensit_ieee64_sum_frst_sp_2d

! second local part of order insensitive summation - dp-variant 
! (sp summation uses this, too; cast is done on invocator side)
! convert back to float-type from scaled integers
  FUNCTION order_insensit_ieee64_sum_scnd_2d(isum, abs_max) &
    & RESULT(rsum)
    INTEGER(KIND=i8), INTENT(IN) :: isum(2)
    REAL(KIND=dp), INTENT(IN) :: abs_max
    REAL(KIND=dp) :: r_fact, rsum
    INTEGER(KIND=i8) :: ival1, ival2
    INTEGER :: iexp
    REAL(dp), PARAMETER :: two_30 = 1073741824._dp
    REAL(dp), PARAMETER :: r_two_30 = 1._dp / two_30
#if defined (__SX__) || defined (__PGI)
    INTEGER(i8) :: mask30
    DATA mask30 / z'000000003fffffff' /
#else
    INTEGER(i8), PARAMETER :: mask30 = INT(z'000000003fffffff',i8)
#endif

    iexp = EXPONENT(abs_max)
    IF (iexp < -980) THEN
      rsum = 0._dp
      RETURN
    END IF
    r_fact = SCALE(1._dp,iexp-30)
    IF (isum(1) >= 0_i8) THEN
      ival1 = ISHFT(isum(1),-30)
      ival2 = IAND (isum(1),mask30)
      rsum = (REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact
    ELSE
      ival1 = ISHFT(ABS(isum(1)),-30)
      ival2 = IAND (ABS(isum(1)),mask30)
      rsum = - ((REAL(ival1,dp)*r_fact)*two_30 + REAL(ival2,dp)*r_fact)
    END IF
    IF (isum(2) >= 0_i8) THEN
      ival1 = ISHFT(isum(2),-30)
      ival2 = IAND (isum(2),mask30)
      rsum = rsum + (REAL(ival1,dp)*r_fact) + (REAL(ival2,dp)*r_fact)*r_two_30
    ELSE
      ival1 = ISHFT(ABS(isum(2)),-30)
      ival2 = IAND (ABS(isum(2)),mask30)
      rsum = rsum - (REAL(ival1,dp)*r_fact) - (REAL(ival2,dp)*r_fact)*r_two_30
    END IF
  END FUNCTION order_insensit_ieee64_sum_scnd_2d

END MODULE mo_ocean_solve_transfer
