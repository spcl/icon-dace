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
! contains extension to solver backend type: BiCG stabilized

MODULE mo_ocean_solve_bicgStab

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend
  USE mo_fortran_tools, ONLY: set_acc_host_or_device
 
  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_ocean_solve_bicgStab
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_bicgStab'

  TYPE, EXTENDS(t_ocean_solve_backend) :: t_ocean_solve_bicgStab
    PRIVATE
! arrays only used by BCGS
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: r0_wp, r_wp, v_wp, &
      & p_wp, ta1_wp, s_wp, ta2_wp
#ifdef __INTEL_COMPILER
!DIR$ ATTRIBUTES ALIGN : 64 :: r0_wp, r_wp, v_wp, p_wp, ta1_wp, s_wp, ta2_wp
#endif
! interfaces
  CONTAINS
    PROCEDURE :: doit_wp => ocean_solve_bicgStab_cal_wp ! override deferred
    PROCEDURE :: doit_sp => ocean_solve_bicgStab_cal_sp ! override deferred
    PROCEDURE, PRIVATE :: recover_arrays_wp => ocean_solve_bicgStab_recover_arrays_wp
    GENERIC, PRIVATE :: recover_arrays => recover_arrays_wp
  END TYPE t_ocean_solve_bicgStab

CONTAINS

! get solver arrays (alloc them, if not done so, yet) - wp-variant  
  SUBROUTINE ocean_solve_bicgStab_recover_arrays_wp(this, x, b, &
      & r0, r, v, p, ta1, s, ta2)
    CLASS(t_ocean_solve_bicgStab), INTENT(INOUT), TARGET :: this
    REAL(KIND=wp), INTENT(OUT), POINTER, DIMENSION(:,:) :: &
      & x, b, r0, r, v, p, ta1, s, ta2

    IF (.NOT.ALLOCATED(this%r_wp)) THEN
      ALLOCATE(this%r_wp(this%trans%nidx, this%trans%nblk), &
        & this%r0_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%v_wp(this%trans%nidx, this%trans%nblk), &
        & this%p_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%ta1_wp(this%trans%nidx, this%trans%nblk), &
        & this%s_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%ta2_wp(this%trans%nidx, this%trans%nblk))
      this%s_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
      this%r0_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
      this%p_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
    END IF
    x => this%x_wp
    b => this%b_wp
    r0 => this%r0_wp
    r => this%r_wp
    v => this%v_wp
    p => this%p_wp
    ta1 => this%ta1_wp
    s => this%s_wp
    ta2 => this%ta2_wp
  END SUBROUTINE ocean_solve_bicgStab_recover_arrays_wp

! actual BiCG-Stab solve (vanilla) - wp-variant
  SUBROUTINE ocean_solve_bicgStab_cal_wp(this, lacc)
    CLASS(t_ocean_solve_bicgStab), INTENT(INOUT) :: this
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    REAL(KIND=wp) :: alpha, beta, omega, tol, tol2
    REAL(KIND=wp) :: rh_glob, rh_glob_o, r0v_glob, ts_glob, tt_glob, rn
    INTEGER :: nidx_e, nblk, iblk, k, m, k_final
    REAL(KIND=wp), POINTER, DIMENSION(:,:), CONTIGUOUS :: &
      & x, b, r0, r, v, p, ta1, s, ta2
    LOGICAL :: done, lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish(this_mod_name, 'OpenACC version currently not tested/validated')
#endif

! retrieve extends of vector to solve
    nblk = this%trans%nblk
    nidx_e = this%trans%nidx_e
    m = this%par%m
    k_final = -1
! retrieve arrays
    CALL this%recover_arrays(x, b, r0, r, v, p, ta1, s, ta2)
    b(nidx_e+1:, nblk) = 0._wp
! compute initial residual and auxiliary vectors
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, ta1)
    ta1(nidx_e+1:, nblk) = 0._wp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, nblk
      r0(:, iblk) = b(:, iblk) - ta1(:, iblk)
      r(:, iblk) = r0(:, iblk)
      p(:, iblk) = r0(:, iblk)
      ta2(:, iblk) = r0(:, iblk) * r0(:, iblk)
      v(:, iblk) = 0._wp
    END DO
!ICON_OMP END PARALLEL DO
    CALL this%trans%global_sum(ta2, rn)
! tolerance
    tol = this%abs_tol_wp
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * SQRT(rn)
      this%abs_tol_wp = tol
    END IF
    tol2 = tol * tol
    done = .false.
    alpha = 1._wp
    rh_glob_o = 1._wp
    omega = 1._wp
! enter BiCG (stab) iteration
    DO k = 1, m - 1
! check if done
      IF (done) CYCLE
! BiCG part
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        ta1(:, iblk) = r0(:, iblk) * r(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      ta1(nidx_e+1:, nblk) = 0._wp
! compute residual norm
      ta2(nidx_e+1:, nblk) = 0._wp
      CALL this%trans%global_sum(ta1, rh_glob, ta2, rn)
! check if already reached desired tolerance
      IF (rn .LE. tol2) THEN
        done = .true.
        k_final = k
        CYCLE
      END IF
! correct search direction (in direction of modified/preconditioned gradient) / update search vector
      beta = (rh_glob / rh_glob_o) * alpha/omega
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        p(:, iblk) = r(:, iblk) + &
          & beta * (p(:, iblk) - omega * v(:,iblk))
      END DO
!ICON_OMP END PARALLEL DO
! bounds exchange search vector + new lhs
      CALL this%trans%sync(p)
      CALL this%lhs%apply(p, v)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        ta1(:, iblk) = r0(:, iblk) * v(:,iblk)
      END DO
!ICON_OMP END PARALLEL DO
! bounds exchange search vector + new lhs
      ta1(nidx_e+1:, nblk) = 0._wp
      CALL this%trans%global_sum(ta1, r0v_glob)
      alpha = rh_glob/r0v_glob
! update BiCG guess
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO  iblk = 1, nblk
        x(:,iblk) = x(:,iblk) + alpha * p(:,iblk)
        s(:,iblk) = r(:, iblk) - alpha * v(:,iblk)
      END DO
!ICON_OMP END PARALLEL DO
! BiCG is done
! now do a GMRES(1), the stabilizer
! bounds exchange search vector + new lhs
      CALL this%trans%sync(s)
      CALL this%lhs%apply(s, r)
! perform a single arnoldi iteration
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO  iblk = 1, nblk
        ta1(:,iblk) = r(:,iblk) * s(:,iblk)
        ta2(:,iblk) = r(:,iblk) * r(:,iblk)
      END DO
!ICON_OMP END PARALLEL DO
      ta1(nidx_e+1:, nblk) = 0._wp
      ta2(nidx_e+1:, nblk) = 0._wp
      CALL this%trans%global_sum(ta1, ts_glob, ta2, tt_glob)
      omega = ts_glob/tt_glob
! update GMRES guess
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO  iblk = 1, nblk
        x(:,iblk) = x(:,iblk) + omega * s(:,iblk)
        r(:,iblk) = s(:,iblk) - omega * r(:,iblk)
      END DO
!ICON_OMP END PARALLEL DO
! bounds exchange search vector + new lhs
      CALL this%trans%sync(x)
      CALL this%lhs%apply(x, ta1)
! prepare residual norm for next step
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO  iblk = 1, nblk
        ta2(:,iblk) = b(:,iblk) - ta1(:, iblk)
        ta2(:,iblk) = ta2(:,iblk) * ta2(:,iblk)
      END DO
!ICON_OMP END PARALLEL DO
      rh_glob_o = rh_glob
    END DO
    this%niter_cal(1) = k_final
    this%res_wp(1) = SQRT(rn)
  END SUBROUTINE ocean_solve_bicgStab_cal_wp

! we should not get here...
  SUBROUTINE ocean_solve_bicgStab_cal_sp(this)
    CLASS(t_ocean_solve_bicgStab), INTENT(INOUT) :: this

    CALL finish(TRIM(this_mod_name)//":ocean_solve_bicgStab_cal_sp()", &
      & "single precision BiCG-Stab solver not implemented")
  END SUBROUTINE ocean_solve_bicgStab_cal_sp

END MODULE mo_ocean_solve_bicgStab
