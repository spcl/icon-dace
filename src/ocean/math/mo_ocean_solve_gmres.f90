! contains extension to solver backend type: GMRES
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





MODULE mo_ocean_solve_gmres
  USE mo_kind, ONLY: sp, wp
  USE mo_exception, ONLY: finish
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend
  USE mo_fortran_tools, ONLY: set_acc_host_or_device
 
  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_ocean_solve_gmres
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_gmres'

  TYPE, EXTENDS(t_ocean_solve_backend) :: t_ocean_solve_gmres
    PRIVATE
! arrays only used by GMRES
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:,:,:) :: v_sp
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:,:) :: w_sp, z_sp, h_sp
    REAL(KIND=sp), ALLOCATABLE, DIMENSION(:) :: s_sp, c_sp, res_t_sp
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:,:) :: v_wp
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:,:) :: w_wp, z_wp, h_wp
    REAL(KIND=wp), ALLOCATABLE, DIMENSION(:) :: s_wp, c_wp, res_t_wp




! interfaces
  CONTAINS
    PROCEDURE :: doit_wp => ocean_solve_gmres_cal_wp ! override deferred
    PROCEDURE :: doit_sp => ocean_solve_gmres_cal_sp ! override deferred
    PROCEDURE, PRIVATE :: recover_arrays_wp => ocean_solve_gmres_recover_arrays_wp
    PROCEDURE, PRIVATE :: recover_arrays_sp => ocean_solve_gmres_recover_arrays_sp
    GENERIC, PRIVATE :: recover_arrays => recover_arrays_wp, recover_arrays_sp
  END TYPE t_ocean_solve_gmres

CONTAINS

! get solver arrays (alloc them, if not done so, yet) - wp-variant
  SUBROUTINE ocean_solve_gmres_recover_arrays_wp(this, v, x, b, w, z, &
    & h, s, c, res)
    IMPLICIT NONE
    CLASS(t_ocean_solve_gmres), TARGET, INTENT(INOUT) :: this
    REAL(KIND=wp), INTENT(OUT), POINTER :: v(:,:,:), x(:,:), &
      & b(:,:), w(:,:), z(:,:), h(:,:), s(:), c(:), res(:)

    IF (.NOT.ALLOCATED(this%z_wp)) THEN
      ALLOCATE(this%s_wp(this%par%m), this%c_wp(this%par%m), &
        this%h_wp(this%par%m,this%par%m), this%z_wp(this%trans%nidx,this%trans%nblk), &
        this%w_wp(this%trans%nidx,this%trans%nblk_a), this%res_t_wp(this%par%m), &
        this%v_wp(this%trans%nidx,this%trans%nblk_a,this%par%m))
      this%w_wp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
      this%v_wp(:, this%trans%nblk:this%trans%nblk_a, :) = 0._wp
    END IF
    v => this%v_wp
    x => this%x_wp
    b => this%b_wp
    w => this%w_wp
    z => this%z_wp
    h => this%h_wp
    s => this%s_wp
    c => this%c_wp
    res => this%res_t_wp
  END SUBROUTINE ocean_solve_gmres_recover_arrays_wp

! actual GMRES-R(n) implementation
  SUBROUTINE ocean_solve_gmres_cal_wp(this, lacc)
    CLASS(t_ocean_solve_gmres), INTENT(INOUT) :: this
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    REAL(wp) :: tol, ci, h_aux
    INTEGER :: jb, nblk, nidx_e, i, k, i_final
    REAL(KIND=wp), POINTER, CONTIGUOUS :: v(:,:,:), x(:,:), b(:,:), &
      & w(:,:), z(:,:), h(:,:), s(:), c(:), res(:), vi(:,:)
    LOGICAL :: done, lzacc

    CALL set_acc_host_or_device(lzacc, lacc)





! set pointers to internal solver-arrays
    CALL this%recover_arrays(v, x, b, w, z, h, s, c, res)
! local domain size
    nblk = this%trans%nblk
    nidx_e = this%trans%nidx_e
! apply (sync+) lhs for initial residuum
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, w)
! compute initial residuum
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) 
     DO jb = 1, nblk
       w(:,jb) = b(:,jb) - w(:,jb)
       z(:,jb) = w(:,jb) * w(:,jb)
    ENDDO
!ICON_OMP END PARALLEL DO
    w(nidx_e + 1:, nblk) = 0._wp
    z(nidx_e + 1:, nblk) = 0._wp
! reduce to residual norm
    CALL this%trans%global_sum(z, h_aux)
    h_aux = SQRT(h_aux)
    res(1) = h_aux
! check, if something to do
    IF (res(1) .EQ. 0._wp) THEN ! already done
      this%niter_cal(1) = 0
      RETURN
    ENDIF
! tolerance
    tol = this%abs_tol_wp
    IF (.NOT.this%par%use_atol) THEN
      tol = this%par%tol * res(1)
      this%abs_tol_wp = tol
    END IF
    done = (ABS(res(1)) < tol)
! enter Arnoldi loop
    DO i = 1, this%par%m - 1
      IF (done) CYCLE
! compute i-th Krylov-vector
      ci = 1.0_wp/h_aux
      vi => v(:,:,i)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO jb = 1, nblk
        vi(:, jb) = w(:, jb) * ci
      END DO
!ICON_OMP END PARALLEL DO
! apply (sync+) lhs
      CALL this%trans%sync(vi)
      CALL this%lhs%apply(vi, w)
      vi(nidx_e + 1:, nblk) = 0._wp
! perform Gram-Schmid orthogonalization
      DO k = 1, i
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO jb = 1, nblk
          z(:,jb) = w(:,jb)*v(:,jb,k)
        END DO
!ICON_OMP END PARALLEL DO
        z(nidx_e + 1:, nblk) = 0._wp
        CALL this%trans%global_sum(z, h(k,i))
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO jb = 1, nblk
          w(:,jb) = w(:,jb) - h(k,i)*v(:,jb,k)
        END DO
!ICON_OMP END PARALLEL DO
      END DO
! next element of Hessenberg-matrix
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO jb = 1, nblk
        z(:,jb) = w(:,jb) * w(:,jb)
      ENDDO
!ICON_OMP END PARALLEL DO
      z(nidx_e + 1:, nblk) = 0._wp
      CALL this%trans%global_sum(z, h(i+1,i))
      h(i+1,i) = SQRT(h(i+1,i))
! apply Givens rotation (1st part)
      h_aux = h(i+1,i)
      DO k = 1, i-1
        ci = h(k  ,i)
        h(k  ,i) =  c(k)*ci + s(k)*h(k+1,i)
        h(k+1,i) = -s(k)*ci + c(k)*h(k+1,i)
      END DO
! compute next (i-th) rotation
      ci = SQRT(h(i,i)*h(i,i) + h_aux*h_aux)
      c(i) = h(i,i) / ci
      s(i) = h_aux / ci
! apply Givens rotation (2nd part)
      h(i,i) = c(i)*h(i,i) + s(i)*h_aux
! new residual norm
      res(i+1) = -s(i)*res(i)
      res(i  ) =  c(i)*res(i)
! check if done
      done = ABS(res(i+1)) < tol
      IF (done) i_final = i + 1
    END DO
! if solve did not converge to desired accuracy (yet) 
    this%niter_cal(1) = MERGE(i_final, -1, done)
    IF (.NOT.done) i_final = this%par%m
! compute Krylov-expansion coeffs
    DO i = i_final - 1, 1, -1
      ci = res(i)
      DO k = i_final - 1, i + 1, -1
        ci = ci - h(i,k)*c(k)
      END DO
      c(i) = ci / h(i,i)
    END DO
! perform Krylov-expansion
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i)
    DO jb = 1, nblk
      DO i = 1, i_final - 1
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      END DO
    END DO
!ICON_OMP END PARALLEL DO
! the returned resdidual norm....
    this%res_wp(1) = ABS(res(i_final))
  END SUBROUTINE ocean_solve_gmres_cal_wp

! same as ocean_solve_gmres_recover_arrays_wp but for single-precision
  SUBROUTINE ocean_solve_gmres_recover_arrays_sp(this, v, x, b, w, z, &
    & h, s, c, res)
    IMPLICIT NONE
    CLASS(t_ocean_solve_gmres), TARGET, INTENT(INOUT) :: this
    REAL(KIND=sp), INTENT(OUT), POINTER :: v(:,:,:), x(:,:), &
      & b(:,:), w(:,:), z(:,:), h(:,:), s(:), c(:), res(:)

    IF (.NOT.ALLOCATED(this%z_sp)) THEN
      ALLOCATE(this%s_sp(this%par_sp%m), this%c_sp(this%par_sp%m), &
        this%res_t_sp(this%par_sp%m), this%h_sp(this%par_sp%m,this%par_sp%m), &
        this%z_sp(this%trans%nidx,this%trans%nblk), &
        this%w_sp(this%trans%nidx,this%trans%nblk_a), &
        this%v_sp(this%trans%nidx,this%trans%nblk_a,this%par_sp%m))
      this%w_sp(:, this%trans%nblk:this%trans%nblk_a) = 0._wp
      this%v_sp(:, this%trans%nblk:this%trans%nblk_a, :) = 0._wp
    END IF
    v => this%v_sp
    x => this%x_sp
    b => this%b_sp
    w => this%w_sp
    z => this%z_sp
    h => this%h_sp
    s => this%s_sp
    c => this%c_sp
    res => this%res_t_sp
  END SUBROUTINE ocean_solve_gmres_recover_arrays_sp

! same as ocean_solve_gmres_cal_wp but for single-precision
  SUBROUTINE ocean_solve_gmres_cal_sp(this)
    CLASS(t_ocean_solve_gmres), INTENT(INOUT) :: this
    REAL(sp) :: tol, ci, h_aux
    INTEGER :: jb, nblk, nidx_e, i, k, i_final
    REAL(KIND=sp), POINTER, CONTIGUOUS :: v(:,:,:), x(:,:), b(:,:), &
      & w(:,:), z(:,:), h(:,:), s(:), c(:), res(:), vi(:,:)
    LOGICAL :: done

    CALL this%recover_arrays(v, x, b, w, z, h, s, c, res)
    nblk = this%trans%nblk
    nidx_e = this%trans%nidx_e
    CALL this%trans%sync(x)
    CALL this%lhs%apply(x, w)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
     DO jb = 1, nblk
       w(:,jb) = b(:,jb) - w(:,jb)
       z(:,jb) = w(:,jb) * w(:,jb)
    ENDDO
!ICON_OMP END PARALLEL DO
    w(nidx_e + 1:, nblk) = 0._sp
    z(nidx_e + 1:, nblk) = 0._sp
    CALL this%trans%global_sum(z, h_aux)
    h_aux = SQRT(h_aux)
    res(1) = h_aux
    IF (res(1) .EQ. 0._sp) THEN ! already done
      this%niter_cal(2) = 0
      RETURN
    ENDIF
    tol = REAL(this%abs_tol_sp, sp)
    IF (.NOT.this%par_sp%use_atol) THEN
      tol = REAL(this%par_sp%tol) * res(1)
      this%abs_tol_sp = REAL(tol, wp)
    END IF
    done = (ABS(res(1)) < tol)
    DO i = 1, this%par_sp%m - 1
      IF (done) CYCLE
      ci = 1.0_sp/h_aux
      vi => v(:,:,i)
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO jb = 1, nblk
        vi(:, jb) = w(:, jb) * ci
      END DO
!ICON_OMP END PARALLEL DO
      CALL this%trans%sync(vi)
      CALL this%lhs%apply(vi, w)
      vi(nidx_e + 1:, nblk) = 0._sp
      DO k = 1, i
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO jb = 1, nblk
          z(:,jb) = w(:,jb)*v(:,jb,k)
        END DO
!ICON_OMP END PARALLEL DO
        z(nidx_e + 1:, nblk) = 0._sp
        CALL this%trans%global_sum(z, h(k,i))
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
        DO jb = 1, nblk
          w(:,jb) = w(:,jb) - h(k,i)*v(:,jb,k)
        END DO
!ICON_OMP END PARALLEL DO
      END DO
!ICON_OMP PARALLEL DO SCHEDULE(STATIC)
      DO jb = 1, nblk
        z(:,jb) = w(:,jb) * w(:,jb)
      ENDDO
!ICON_OMP END PARALLEL DO
      z(nidx_e + 1:, nblk) = 0._sp
      CALL this%trans%global_sum(z, h(i+1,i))
      h(i+1,i) = SQRT(h(i+1,i))
      h_aux = h(i+1,i)
      DO k = 1, i-1
        ci = h(k  ,i)
        h(k  ,i) =  c(k)*ci + s(k)*h(k+1,i)
        h(k+1,i) = -s(k)*ci + c(k)*h(k+1,i)
      END DO
      ci = SQRT(h(i,i)*h(i,i) + h_aux*h_aux)
      c(i) = h(i,i) / ci
      s(i) = h_aux / ci
      h(i,i) = c(i)*h(i,i) + s(i)*h_aux
      res(i+1) = -s(i)*res(i)
      res(i  ) =  c(i)*res(i)
      done = ABS(res(i+1)) < tol
      IF (done) i_final = i + 1
    END DO
    IF (.NOT.done) i_final = this%par_sp%m
    this%niter_cal(2) = MERGE(i_final, -1, done)
    DO i = i_final - 1, 1, -1
      ci = res(i)
      DO k = i_final - 1, i + 1, -1
        ci = ci - h(i,k)*c(k)
      END DO
      c(i) = ci / h(i,i)
    END DO
!ICON_OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i)
    DO jb = 1, nblk
      DO i = 1, i_final - 1
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      END DO
    END DO
!ICON_OMP END PARALLEL DO
    this%res_wp(2) = REAL(ABS(res(i_final)), wp)
    CALL this%trans%sync(x)
  END SUBROUTINE ocean_solve_gmres_cal_sp

END MODULE mo_ocean_solve_gmres
