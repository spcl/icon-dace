! contains extension to solver backend type: CG + Jacobi preconditioner
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





MODULE mo_ocean_solve_cgj

  USE mo_kind, ONLY: wp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend
  USE mo_fortran_tools, ONLY: set_acc_host_or_device
 
  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_ocean_solve_cgj
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_cgj'

  TYPE, EXTENDS(t_ocean_solve_backend) :: t_ocean_solve_cgj
    PRIVATE
! arrays only used by CGJ
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: &
      & z_wp, d_wp, r1_wp, rsq_wp, h_wp



! interfaces
  CONTAINS
    PROCEDURE :: doit_wp => ocean_solve_cgj_cal_wp ! override deferred
    PROCEDURE :: doit_sp => ocean_solve_cgj_cal_sp ! override deferred
    PROCEDURE, PRIVATE :: recover_arrays_wp => ocean_solve_cgj_recover_arrays_wp
    GENERIC, PRIVATE :: recover_arrays => recover_arrays_wp
  END TYPE t_ocean_solve_cgj

CONTAINS

! get solver arrays (alloc them, if not done so, yet) - wp-variant  
SUBROUTINE ocean_solve_cgj_recover_arrays_wp(this, x, b, z, d, r, r2, &
    & h)
    CLASS(t_ocean_solve_cgj), INTENT(INOUT), TARGET :: this
    REAL(KIND=wp), INTENT(OUT), POINTER, DIMENSION(:,:) :: &
      & x, b, z, d, r, r2, h

    IF (.NOT.ALLOCATED(this%z_wp)) THEN
      ALLOCATE(this%z_wp(this%trans%nidx, this%trans%nblk), &
        & this%d_wp(this%trans%nidx, this%trans%nblk_a), &
        & this%r1_wp(this%trans%nidx, this%trans%nblk), &
        & this%rsq_wp(this%trans%nidx, this%trans%nblk), &
        & this%h_wp(this%trans%nidx, this%trans%nblk))
      this%d_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
    END IF
    x => this%x_wp
    b => this%b_wp
    z => this%z_wp
    d => this%d_wp
    r => this%r1_wp
    r2 => this%rsq_wp
    h => this%h_wp
  END SUBROUTINE ocean_solve_cgj_recover_arrays_wp

! actual CG solve utilizing Jacobi preconditioner - wp-variant
  SUBROUTINE ocean_solve_cgj_cal_wp(this, lacc)
    CLASS(t_ocean_solve_cgj), INTENT(INOUT) :: this
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    REAL(KIND=wp) :: alpha, beta, dz_glob, tol, tol2
    REAL(KIND=wp) :: rh_glob, rh_glob_o, rn
    INTEGER :: nidx_a, nidx_e, nblk, iblk, k, m, k_final
    REAL(KIND=wp), POINTER, DIMENSION(:,:), CONTIGUOUS :: &
      & x, b, z, d, r, r2, invaii, h
    LOGICAL :: done, lzacc

    CALL set_acc_host_or_device(lzacc, lacc)





! retrieve extends of vector to solve
    nidx_a = this%trans%nidx
    nblk = this%trans%nblk
    nidx_e = this%trans%nidx_e
    m = this%par%m
    k_final = -1
! retrieve arrays
    CALL this%recover_arrays(x, b, z, d, r, r2, h)
! retrieve preconditioner matrix (i.e. Jij = 1/Aij , if i=j and Jij = 0. otherwise
    CALL this%lhs%get_invaii(invaii)
    b(nidx_e+1:, nblk) = 0._wp
! compute initial residual and auxiliary vectors
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, z)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
    DO iblk = 1, nblk
      r(:, iblk) = b(:, iblk) - z(:, iblk)
      r2(:, iblk) = r(:, iblk) * r(:, iblk)
      h(:, iblk) = r(:, iblk) * invaii(:, iblk)
      z(:, iblk) = r(:, iblk) * h(:, iblk)
      d(:, iblk) = h(:, iblk)
    END DO
!ICON_OMP END PARALLEL DO
    CALL this%trans%global_sum(r2, rn, z, rh_glob)
! tolerance
    tol = this%abs_tol_wp
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * SQRT(rn)
      this%abs_tol_wp = tol
    END IF
    tol2 = tol * tol
    done = .false.
! enter CG iteration
    DO k = 1, m - 1
! check if done
      IF (done) CYCLE
! compute residual norm
      r2(nidx_e+1:nidx_a, nblk) = 0._wp
      z(nidx_e+1:, nblk) = 0._wp
      IF (k .GT. 1) CALL this%trans%global_sum(r2, rn, z, rh_glob)
! check if already reached desired tolerance
      IF (rn .LE. tol2) THEN
        done = .true.
        k_final = k
        CYCLE
      END IF
! correct search direction (in direction of modified/preconditioned gradient) / update search vector
      IF (k .GT. 1) THEN
        beta = rh_glob / rh_glob_o
        d(nidx_e+1:, nblk) = 0._wp
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO iblk = 1, nblk
          d(:, iblk) = h(:, iblk) + beta * d(:, iblk)
        END DO
!ICON_OMP END PARALLEL DO
      END IF
! bounds exchange search vector + new lhs
      CALL this%trans%sync(d)
      CALL this%lhs%apply(d, z)
      d(nidx_e+1:, nblk) = 0._wp
! compute extrapolated location of minimum in direction of d
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        r2(:, iblk) = d(:, iblk) * z(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      r2(nidx_e+1:, nblk) = 0._wp
      CALL this%trans%global_sum(r2, dz_glob)
      alpha = rh_glob / dz_glob
! update guess and residuum, apply preconditioner
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO iblk = 1, nblk
        x(:, iblk) = x(:, iblk) + alpha * d(:, iblk)
        r(:, iblk) = r(:, iblk) - alpha * z(:, iblk)
        r2(:, iblk) = r(:, iblk) * r(:, iblk)
        h(:, iblk) = r(:, iblk) * invaii(:, iblk)
        z(:, iblk) = r(:, iblk) * h(:, iblk)
      END DO
!ICON_OMP END PARALLEL DO
      rh_glob_o = rh_glob
    END DO
    this%niter_cal(1) = k_final
    this%res_wp(1) = SQRT(rn)
    CALL this%trans%sync(x)
  END SUBROUTINE ocean_solve_cgj_cal_wp

! we should not get here...
  SUBROUTINE ocean_solve_cgj_cal_sp(this)
    CLASS(t_ocean_solve_cgj), INTENT(INOUT) :: this

    CALL finish(TRIM(this_mod_name)//":ocean_solve_cgj_cal_sp()", &
      & "single precision CG solver with Jacobi-preconditioner not implemented")
  END SUBROUTINE ocean_solve_cgj_cal_sp

END MODULE mo_ocean_solve_cgj
