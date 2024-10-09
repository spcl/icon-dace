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
! contains extension to solver backend type: MINRES

MODULE mo_ocean_solve_mres

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend
  USE mo_fortran_tools, ONLY: set_acc_host_or_device
 
  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_ocean_solve_mres
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_mres'

  TYPE, EXTENDS(t_ocean_solve_backend) :: t_ocean_solve_mres
    PRIVATE
! arrays only used by MINRES
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: p0_wp, p1_wp, p2_wp, s0_wp, &
      & s1_wp, s2_wp, r_wp, ss_wp, rs_wp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: p0_wp, p1_wp, p2_wp, s0_wp, s1_wp, s2_wp, r_wp
!DIR$ ATTRIBUTES ALIGN : 64 :: ss_wp, rs_wp
#endif
! interfaces
  CONTAINS
    PROCEDURE :: doit_wp => ocean_solve_mres_cal_wp ! override deferred
    PROCEDURE :: doit_sp => ocean_solve_mres_cal_sp ! override deferred
    PROCEDURE, PRIVATE :: recover_arrays_wp => ocean_solve_mres_recover_arrays_wp
    GENERIC, PRIVATE :: recover_arrays => recover_arrays_wp
  END TYPE t_ocean_solve_mres

CONTAINS

! get solver arrays (alloc them, if not done so, yet) - wp-variant
  SUBROUTINE ocean_solve_mres_recover_arrays_wp(this, x, b, r, p0, &
    & p1, p2, s0, s1, s2, rs, ss)
    CLASS(t_ocean_solve_mres), INTENT(INOUT), TARGET :: this
    REAL(KIND=wp), INTENT(INOUT), POINTER, DIMENSION(:,:) :: &
      & x, b, r, p0, p1, p2, s0, s1, s2, rs, ss

    IF (.NOT.ALLOCATED(this%r_wp)) THEN
      ALLOCATE(this%r_wp(this%trans%nidx, this%trans%nblk), &
        & this%p0_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%p1_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%p2_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%s0_wp(this%trans%nidx, this%trans%nblk), &
        & this%s1_wp(this%trans%nidx, this%trans%nblk), &
        & this%s2_wp(this%trans%nidx, this%trans%nblk), &
        & this%rs_wp(this%trans%nidx, this%trans%nblk), &
        & this%ss_wp(this%trans%nidx, this%trans%nblk))
      this%p0_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
      this%p1_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
      this%p2_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
    END IF
    x => this%x_wp
    b => this%b_wp
    r => this%r_wp
    p0 => this%p0_wp
    p1 => this%p1_wp
    p2 => this%p2_wp
    s0 => this%s0_wp
    s1 => this%s1_wp
    s2 => this%s2_wp
    ss => this%ss_wp
    rs => this%rs_wp
  END SUBROUTINE ocean_solve_mres_recover_arrays_wp

! actual MINRES solve (vanilla) - wp-variant
SUBROUTINE ocean_solve_mres_cal_wp(this, lacc)
    CLASS(t_ocean_solve_mres), INTENT(INOUT) :: this
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    REAL(KIND=wp) :: alpha, beta1, beta2, tol, tol2, rn, rsg, ssg, ss0g, ssg_o
    INTEGER :: nidx_a, nidx_e, nblk, iblk, k, m, k_final
    REAL(KIND=wp), POINTER, DIMENSION(:,:), CONTIGUOUS :: &
      & x, b, tmp, r, p0, p1, p2, s0, s1, s2, rs, ss
    LOGICAL :: done, lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish(this_mod_name, 'OpenACC version currently not tested/validated')
#endif

! retrieve extends of vector to solve
    nidx_a = this%trans%nidx
    nblk = this%trans%nblk
    nidx_e = this%trans%nidx_e
    m = this%par%m
    k_final = -1
! retrieve arrays
    CALL this%recover_arrays(x, b, r, p0, p1, p2, s0, s1, s2, rs, ss)
#ifdef __PGI
    ALLOCATE(tmp(this%trans%nidx, this%trans%nblk))
#endif
    b(nidx_e+1:, nblk) = 0._wp
! compute initial residual and auxiliary vectors
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, ss)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, nblk
      r(:, iblk) = b(:, iblk) - ss(:, iblk)
      p0(:, iblk) = r(:, iblk)
      p1(:, iblk) = p0(:, iblk)
      rs(:, iblk) = r(:, iblk) * r(:, iblk)
    END DO
!ICON_OMP END PARALLEL DO
    rs(nidx_e+1:, nblk) = 0._wp
    CALL this%trans%global_sum(rs, rn)
    CALL this%trans%sync(p0)
    CALL this%lhs%apply(p0, s0)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, nblk
      s1(:, iblk) = s0(:, iblk)
      ss(:, iblk) = s0(:, iblk) * s0(:, iblk)
    END DO
!ICON_OMP END PARALLEL DO
    tol = this%abs_tol_wp
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * SQRT(rn)
      this%abs_tol_wp = tol
    END IF
    tol2 = tol * tol
    done = .false.
! enter MINRES iteration
    DO k = 1, m
! check if done
      IF (done) CYCLE
#ifdef __PGI
! pgi does not like pointer juggling !?!
      DO iblk = 1, nblk
        tmp(:, iblk) = p2(:, iblk)
        p2(:, iblk) = p1(:, iblk)
        p1(:, iblk) = p0(:, iblk)
        p0(:, iblk) = tmp(:, iblk)
        tmp(:, iblk) = s2(:, iblk)
        s2(:, iblk) = s1(:, iblk)
        s1(:, iblk) = s0(:, iblk)
        s0(:, iblk) = tmp(:, iblk)
      END DO
#else
      tmp => p2
      p2 => p1
      p1 => p0
      p0 => tmp
      tmp => s2
      s2 => s1
      s1 => s0
      s0 => tmp
#endif
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        rs(:, iblk) = s1(:, iblk) * r(:, iblk)
        ss(:, iblk) = s1(:, iblk) * s1(:, iblk)
        p0(:, iblk) = s1(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      rs(nidx_e+1:, nblk) = 0._wp
      ss(nidx_e+1:, nblk) = 0._wp
      CALL this%trans%global_sum(rs, rsg, ss, ssg)
      alpha = rsg / ssg
      CALL this%trans%sync(p0)
      CALL this%lhs%apply(p0, s0)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        x(:, iblk) = x(:, iblk) + alpha * p1(:, iblk)
        r(:, iblk) = r(:, iblk) - alpha * s1(:, iblk)
        rs(:, iblk) = r(:, iblk) * r(:, iblk)
        ss(:, iblk) = s1(:, iblk) * s0(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
! compute residual norm
      rs(nidx_e+1:nidx_a, nblk) = 0._wp
      ss(nidx_e+1:, nblk) = 0._wp
      CALL this%trans%global_sum(rs, rn, ss, ss0g)
! check if already reached desired tolerance
      IF (rn .LE. tol2 .OR. ABS(alpha) .GT. 1.e100_wp) THEN
        done = .true.
        k_final = k
        CYCLE
      END IF
      beta1 = ss0g / ssg
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        p0(:, iblk) = p0(:, iblk) - beta1 * p1(:, iblk)
        s0(:, iblk) = s0(:, iblk) - beta1 * s1(:, iblk)
        ss(:, iblk) = s2(:, iblk) * s0(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      IF (k .GT. 1) THEN
        ss(nidx_e+1:, nblk) = 0._wp
        CALL this%trans%global_sum(ss, ss0g)
        beta2 = ss0g / ssg_o
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO iblk = 1, nblk
          p0(:, iblk) = p0(:, iblk) - beta2 * p2(:, iblk)
          s0(:, iblk) = s0(:, iblk) - beta2 * s2(:, iblk)
        END DO
!ICON_OMP END PARALLEL DO
      END IF
      ssg_o = ssg
    END DO
    this%niter_cal(1) = k_final
    this%res_wp(1) = SQRT(rn)
#ifdef __PGI
    DEALLOCATE(tmp)
#endif
  END SUBROUTINE ocean_solve_mres_cal_wp

! we should not get here...
  SUBROUTINE ocean_solve_mres_cal_sp(this)
    CLASS(t_ocean_solve_mres), INTENT(INOUT) :: this

    CALL finish(TRIM(this_mod_name)//":ocean_solve_mres_cal_sp()", &
      & "single precision minres solver not implemented")
  END SUBROUTINE ocean_solve_mres_cal_sp

END MODULE mo_ocean_solve_mres
