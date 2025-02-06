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

! Contains extension to solver backend type: GMRES (legacy implementation)

!#include "omp_definitions.inc"

MODULE mo_ocean_solve_legacy_gmres
  USE mo_kind, ONLY: wp
  USE mo_ocean_solve_lhs, ONLY: t_lhs
  USE mo_ocean_solve_transfer, ONLY: t_transfer
  USE mo_ocean_solve_backend, ONLY: t_ocean_solve_backend
  USE mo_mpi, ONLY: p_sum
  USE mo_exception, ONLY: finish
  USE mo_timer, ONLY: timer_start, timer_stop
  USE mo_run_config, ONLY: ltimer
  USE mo_fortran_tools, ONLY: set_acc_host_or_device
 
  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: t_ocean_solve_legacy_gmres
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'mo_ocean_solve_legacy_gmres'

  TYPE, EXTENDS(t_ocean_solve_backend) :: t_ocean_solve_legacy_gmres
    PRIVATE
! interfaces
  CONTAINS
    PROCEDURE :: doit_wp => ocean_solve_legacy_gmres_cal_wp ! override deferred
    PROCEDURE :: doit_sp => ocean_solve_legacy_gmres_cal_sp ! override deferred
  END TYPE t_ocean_solve_legacy_gmres

CONTAINS

! actual GMRES solve (vanilla)
  SUBROUTINE ocean_solve_legacy_gmres_cal_wp(this, lacc)
    CLASS(t_ocean_solve_legacy_gmres), INTENT(INOUT) :: this
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    LOGICAL :: maxiterex ! is reconstructed later
    LOGICAL :: lzacc
    INTEGER :: niter

    CALL set_acc_host_or_device(lzacc, lacc)

#ifdef _OPENACC
    IF (lzacc) CALL finish(this_mod_name, 'OpenACC version currently not tested/validated')
#endif

    CALL ocean_restart_gmres(this%x_wp, this%lhs, this%b_wp, &
      & this%par%tol, this%par%use_atol, this%par%m, maxiterex, &
      & niter, this%res_wp, this%trans)
    this%niter_cal(1) = MERGE(niter, -1, .NOT.maxiterex)
  END SUBROUTINE ocean_solve_legacy_gmres_cal_wp

! same as ocean_solve_legacy_gmres_cal_wp but for single-precision
  SUBROUTINE ocean_solve_legacy_gmres_cal_sp(this)
    CLASS(t_ocean_solve_legacy_gmres), INTENT(INOUT) :: this
    CHARACTER(len=*), PARAMETER :: routine='ocean_solve_legacy_gmres'
 
    CALL finish(routine, "not implemented!")
  END SUBROUTINE ocean_solve_legacy_gmres_cal_sp

! the legacy backend

  !-------------------------------------------------------------------------
  !>
  !!
  !! Based on gmres_oce_old
  !-------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE ocean_restart_gmres( x,lhs,        &
!    & h_e, old_h, patch_3D,      &
!    & p_op_coeff, &
    & b,                            &
    & tolerance,  use_absolute_tolerance,  &
    & m, maxiterex,niter,res,                   &
!    & , preconditioner ! no preconditioner here
    & trans &
    & )
    !
    ! !DESCRIPTION
    !  GMRES solver for linear systems. Notice that, in principle, r and b
    !  can be of any TYPE for which the following operations are defined:
    !   +, -, dot_product
    !
    REAL(wp), INTENT(inout) :: x(:,:)
    TYPE(t_lhs), INTENT(INOUT) :: lhs
!    REAL(wp), INTENT(in) :: h_e(:,:)
!    REAL(wp), INTENT(in) :: old_h(:,:)
    ! patch info needed for calculating lhs
    !TYPE(t_patch), INTENT(IN) :: patch_2D
!    TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
    ! index defining the "active" region of the arrays
    ! parameter used in calculating the lhs
    ! REAL(wp), INTENT(in) :: coeff
!    TYPE(t_operator_coeff), INTENT(in):: p_op_coeff
    ! right-hand side: same shape as x
    REAL(wp), INTENT(inout) :: b(:,:) ! same size as x
    REAL(wp), INTENT(inout) :: tolerance
    LOGICAL, INTENT(in) :: use_absolute_tolerance ! (relative or absolute) tolerance
    ! .false. for relative relative_tolerance
    INTEGER,  INTENT(in) :: m         ! maximum number of iterations
    LOGICAL,  INTENT(out) :: maxiterex ! true if reached m iterations
    INTEGER,  INTENT(out) :: niter    ! number of iterations (defined
    ! as the number of evaluation of
    ! the lhs which, due to the
    ! initial evaluation of the
    ! residual, is equal to the
    ! number of Arnoldi iterations +1)
    ! norms of the residual (convergence history); an argument of
    ! dimension at least m is required
    REAL(wp), INTENT(inout) :: res(1) ! (m)
    CLASS(t_transfer), POINTER, INTENT(INOUT) :: trans

!    INTERFACE   ! left-hand-side: A*x
!      FUNCTION lhs(x, patch_3D, h_e, p_op_coeff) result(ax)
!        USE mo_kind, ONLY: wp
!        USE mo_model_domain, ONLY: t_patch, t_patch_3d
!        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
!        REAL(wp),    INTENT(inout) :: x(:,:)  ! inout for sync
!!         REAL(wp), INTENT(in) :: old_h(:,:)
!        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
!        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
!        ! REAL(wp),    INTENT(in) :: coeff
!        REAL(wp),    INTENT(in) :: h_e(:,:)
!        TYPE(t_operator_coeff),INTENT(in)  :: p_op_coeff
!        !     REAL(wp), POINTER :: ax( :, : ) ! same as x
!        REAL(wp) :: ax( SIZE(x,1) , SIZE(x,2) ) ! same as x
!      END FUNCTION lhs
!    END INTERFACE

!    INTERFACE   ! preconditioner
!      SUBROUTINE preconditioner(r, patch_3D, p_op_coeff,h_e)
!        USE mo_kind, ONLY: wp
!        USE mo_model_domain, ONLY: t_patch, t_patch_3d
!        USE mo_operator_ocean_coeff_3d, ONLY: t_operator_coeff
!        REAL(wp), INTENT(inout)           :: r(:,:)
!        !TYPE(t_patch), TARGET, INTENT(in) :: patch_2D
!        TYPE(t_patch_3d ),TARGET, INTENT(in)   :: patch_3D
!        TYPE(t_operator_coeff),INTENT(in) :: p_op_coeff
!        REAL(wp),    INTENT(in)           :: h_e(:,:)
!      END SUBROUTINE preconditioner
!    END INTERFACE
!
!    OPTIONAL :: preconditioner

    ! !LOCAL VARIABLES

    LOGICAL :: done
    INTEGER ::     &
      & i,           & ! index for the Arnoldi loop
      & k,           & ! index for the Gram Schmidt orthogonalization
      & nkry           ! dimension of the Krylov space
    REAL(wp) ::    &
      & tol,         & ! effective relative_tolerance (absolute o relative)
      & tol2,        &
      & r(SIZE(x,1),SIZE(x,2)),  & ! residual
      & rn2(m),      &             ! two-norm of the residual
      & v(SIZE(x,1),SIZE(x,2),m),& ! Krylov basis
      & h(m,m),      &             ! Hessemberg matrix
      & hki, hk1i,   &
      & c(m), s(m),  &             ! rotation matrices
      & den, ci

    ! REAL(wp), POINTER :: w(:,:)
    REAL(wp) ::  w(SIZE(x,1),SIZE(x,2)) ! new Krylov vector

    REAL(wp) :: rrn2, h_aux, rh

    INTEGER :: jb, no_of_blocks, pad_nproma, nproma, my_mpi_work_communicator

!    REAL(wp) :: sum_aux(patch_3D%p_patch_2d(1)%cells%in_domain%end_block)
    REAL(wp) :: sum_aux(trans%nblk)
    !   REAL(wp) :: sum_x(0:7), sum_w, sum_v
    ! #else
    !   REAL(wp) :: z(SIZE(x,1),SIZE(x,2)) ! needed for global sums in
    !   p_test_run
    ! #endif

    INTEGER :: mythreadno
!    TYPE(t_patch), POINTER :: patch_2d
    CONTIGUOUS :: b, x

#ifdef _OPENMP
    INTEGER OMP_GET_THREAD_NUM
#endif
    !-------------------------------------------------------------------------
    !  write(0,*) "--------------- gmres --------------------------"

!    patch_2d =>  patch_3D%p_patch_2d(1)
    my_mpi_work_communicator = trans%comm
    ! 0) set module variables and initialize maxiterex
!    no_of_blocks    = patch_2d%cells%in_domain%end_block
!    end_nproma      = patch_2d%cells%in_domain%end_index
!    pad_nproma      = patch_2d%cells%in_domain%end_index + 1
    nproma          = trans%nidx
    no_of_blocks    = trans%nblk
!    end_nproma      = trans%nidx_e
    pad_nproma      = trans%nidx_e + 1

    maxiterex = .FALSE.

    v  = 0.0_wp
    r = 0.0_wp
    b(pad_nproma:nproma, no_of_blocks) = 0.0_wp
 !   w(:, :) = 0.0_wp
 !   b(:, :) = 0.0_wp
 !   x(:, :) = 0.0_wp

    ! 1) compute the preconditioned residual
!    IF (PRESENT(preconditioner)) CALL preconditioner(x(:,:),patch_3D,p_op_coeff,h_e)
!    IF (PRESENT(preconditioner)) CALL preconditioner(b(:,:),patch_3D,p_op_coeff,h_e)


!    w(:, :) = lhs(x(:,:), patch_3D, h_e, p_op_coeff)
    CALL trans%sync(x)
    CALL lhs%apply(x, w, .true.)

    w(pad_nproma:nproma, no_of_blocks) = 0.0_wp
    x(pad_nproma:nproma, no_of_blocks) = 0.0_wp

    !    write(0,*) "-----------------------------------"
    !    sum_x(1) = SUM(x(:,:))
    !    sum_w    = SUM(w(:,:))
    !    write(0,*) "gmres sum x(0:1), w:", sum_x(0:1), sum_w

!    start_timer(timer_gmres,2)

    mythreadno = 0
! !ICON_OMP_PARALLEL PRIVATE(rrn2, myThreadNo)
!    !$ myThreadNo = OMP_GET_THREAD_NUM()

! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      r(1:nproma,jb) = b(1:nproma,jb) - w(1:nproma,jb)
    ENDDO
! !ICON_OMP_END_DO

!     IF (PRESENT(preconditioner)) CALL
!     preconditioner(r(:,:),patch_3D,p_op_coeff,h_e)
    IF (ltimer) CALL timer_start(trans%timer_glob_sum)
! !ICON_OMP_DO PRIVATE(jb) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      sum_aux(jb) = SUM(r(1:nproma,jb) * r(1:nproma,jb))
    ENDDO
    ! sum_aux(no_of_blocks) = SUM(r(1:end_nproma,no_of_blocks) *
    ! r(1:end_nproma,no_of_blocks))
! !ICON_OMP_END_DO

    IF (mythreadno == 0) THEN
      h_aux  = SUM(sum_aux(1:no_of_blocks))
!      start_detail_timer(timer_gmres_p_sum, 5)
      h_aux = p_sum(h_aux, my_mpi_work_communicator)
      IF (ltimer) CALL timer_stop(trans%timer_glob_sum)
!      stop_detail_timer(timer_gmres_p_sum, 5)
      rn2(1) = SQRT(h_aux)
      !       rn2(1) = SQRT(p_sum(h_aux, my_mpi_work_communicator))
      ! !ICON_OMP FLUSH(rn2(1))
    ENDIF
! !ICON_OMP BARRIER


    ! 2) compute the first vector of the Krylov space
    IF (rn2(1) /= 0.0_wp) THEN
      rrn2 = 1.0_wp/rn2(1)
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        v(1:nproma,jb,1) = r(1:nproma,jb) * rrn2
      ENDDO
! !ICON_OMP_END_DO NOWAIT
    ENDIF
! !ICON_OMP_END_PARALLEL

!    CALL dbg_print('1: w', w, method_name, 3, in_subset=patch_2d%cells%owned)
!    CALL dbg_print('1: r', r, method_name, 3, in_subset=patch_2d%cells%owned)
!    write(0,*) "h_aux:", h_aux, " rn2(1):", rn2(1)

    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'r=',
    !      r(pad_nproma+1:nproma, no_of_blocks)
    !      write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'v=',
    !      v(pad_nproma+1:nproma, no_of_blocks, 1)
    !      CALL p_barrier

    IF (rn2(1) == 0.0_wp) THEN ! already done
      !    print*,' gmres: rn2(1)=0.0 '
      ! #slo# 2010-06-15 trying to exit with niter=1 and residual=0.0
      niter = 0
      res(1) = ABS(rn2(1))

!      stop_timer(timer_gmres,2)

      RETURN
    ENDIF

    ! 3) define the relative_tolerance: can be absolute or relative
    IF (use_absolute_tolerance) THEN
      tol = tolerance
    ELSE
      tol = tolerance * rn2(1)
      tolerance = tol
    ENDIF
    tol2 = tol**2

    ! 4) Arnoldi loop
    arnoldi: DO i = 1, m-1
      ! 4.1) compute the next (i.e. i+1) Krylov vector
      !      sum_v = SUM(v(:,:,i))
      !       write(0,*) i, " gmres v before lhs:", sum_v, sum_w
      ! write(*,*)  get_my_global_mpi_id(), 'before lhs v(pad):', pad_nproma,
      ! nproma, k, 'v=', v(pad_nproma:nproma, no_of_blocks, k)
      CALL trans%sync(v(:,:,i))
      call lhs%apply(v(:,:,i), w, .true.)
!      w(:,:) = lhs( v(:,:,i), patch_3D, h_e, p_op_coeff )
      !w(pad_nproma:nproma, no_of_blocks)   = 0.0_wp
      v(pad_nproma:nproma, no_of_blocks,i) = 0.0_wp

      !      sum_v = SUM(v(:,:,i))
      !      sum_w = SUM(w(:,:))
      !      write(0,*) i, " gmres sum v, w:", sum_v, sum_w

      ! 4.2) Gram-Schmidt orthogonalization


!      IF (PRESENT(preconditioner)) CALL preconditioner(w(:,:),patch_3D,p_op_coeff,h_e)

      gs_orth: DO k = 1, i

        IF (ltimer) CALL timer_start(trans%timer_glob_sum)
!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
#ifdef _OPENMP
        myThreadNo = OMP_GET_THREAD_NUM()
#endif
!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          sum_aux(jb) = SUM(w(1:nproma,jb) * v(1:nproma,jb,k))
        ENDDO
!ICON_OMP_END_DO

        ! write(*,*)  get_my_global_mpi_id(), 'v(pad):', pad_nproma, nproma, k,
        ! 'v=', v(pad_nproma:nproma, no_of_blocks, k)
        ! write(*,*)  get_my_global_mpi_id(), 'w(pad):', pad_nproma, nproma,
        ! 'w=', w(pad_nproma:nproma, no_of_blocks)
        ! CALL p_barrier

        IF (mythreadno == 0) THEN
          h_aux = SUM(sum_aux(1:no_of_blocks))
!          start_detail_timer(timer_gmres_p_sum,5)
          h_aux = p_sum(h_aux, my_mpi_work_communicator)
!          stop_detail_timer(timer_gmres_p_sum,5)
          h(k,i) = h_aux
          IF (ltimer) CALL timer_stop(trans%timer_glob_sum)
!ICON_OMP FLUSH(h_aux)
        ENDIF
!ICON_OMP BARRIER
!        CALL dbg_print('anroldi: w', w, method_name, 3,
!        in_subset=patch_2d%cells%owned)
!        CALL dbg_print('anroldi: v', v, method_name, 3,
!        in_subset=patch_2d%cells%owned)
!        write(0,*) "h_aux:", h_aux

!ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          w(1:nproma, jb) = w(1:nproma, jb) - h_aux * v(1:nproma, jb, k)
        ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
      ENDDO gs_orth

      !    write(*,*)  get_my_global_mpi_id(), ':', pad_nproma, nproma, 'w=',
      !    w(pad_nproma+1:nproma, no_of_blocks)
      !    CALL p_barrier
      !    w(pad_nproma:nproma, no_of_blocks) = 0.0_wp
      !    CALL sync_patch_array(SYNC_C, patch_2D, w(:,:) )

! !ICON_OMP_PARALLEL PRIVATE(rh, myThreadNo)
! !$   myThreadNo = OMP_GET_THREAD_NUM()
      ! 4.3) new element for h
      IF (ltimer) CALL timer_start(trans%timer_glob_sum)
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
      DO jb = 1, no_of_blocks
        sum_aux(jb) = SUM(w(1:nproma,jb) * w(1:nproma,jb))
      ENDDO
! !ICON_OMP_END_DO
      IF (mythreadno == 0) THEN
        h_aux = SUM(sum_aux(1:no_of_blocks))
!        start_detail_timer(timer_gmres_p_sum,5)
        h_aux = p_sum(h_aux, my_mpi_work_communicator)
!        stop_detail_timer(timer_gmres_p_sum,5)
        IF (ltimer) CALL timer_stop(trans%timer_glob_sum)
        h_aux = SQRT(h_aux)
        h(i+1,i) = h_aux
        done = h_aux .LT. tol2
      ENDIF
! !ICON_OMP FLUSH(h_aux, done)
! !ICON_OMP BARRIER
      !     write(0,*) i, " gmres  tol2, h_aux", tol2, h_aux


      ! 4.4) if w is independent from v, add v(:,:,:,i+1)
      IF (.NOT. done) THEN
        rh = 1.0_wp/h_aux
! !ICON_OMP_DO ICON_OMP_DEFAULT_SCHEDULE
        DO jb = 1, no_of_blocks
          v(1:nproma, jb, i+1) = w(1:nproma,jb) * rh
        ENDDO
! !ICON_OMP_END_DO NOWAIT
      ENDIF

! !ICON_OMP_END_PARALLEL

      ! 4.5) apply the rotation matrices
      rotation: DO k = 1, i-1
        hki  = h(k  ,i)
        hk1i = h(k+1,i)
        h(k  ,i) =  c(k)*hki + s(k)*hk1i
        h(k+1,i) = -s(k)*hki + c(k)*hk1i
      ENDDO rotation

      ! 4.6) compute the new (i.e. i) rotation
      den = SQRT(h(i,i)**2 + h(i+1,i)**2)
      c(i) = h(i  ,i) / den
      s(i) = h(i+1,i) / den

      ! 4.7) complete applying the rotation matrices
      h(i,i) =  c(i)*h(i,i) + s(i)*h(i+1,i)

      ! 4.8) compute new residual norm
      rn2(i+1) = -s(i)*rn2(i)
      rn2(i  ) =  c(i)*rn2(i)

      ! 4.9) check whether we are done
      ! write(0,*) i, " gmres  relative_tolerance bound, residual: ", tol,
      ! ABS(rn2(i+1))

      IF ( done .OR. (ABS(rn2(i+1)) < tol) ) THEN
       done = .true.
       EXIT arnoldi
      END IF
    ENDDO arnoldi

    ! 5) check whether we are here because we have exceeded m
    !    Notice that if the arnoldi loop has been terminated with the
    !    exit instruction, we have i=m-1, while if the loop has been
    !    terminated by the loop counter we have i=m (i is incremented
    !    after finishing the do loop)
    IF (.NOT. done) THEN
      maxiterex = .TRUE.
      ! only the first m-1 columns h(i,i) have been set in the arnoldi
      ! loop, but we have i=m. We thus decrement i.
      i = i-1
    ENDIF

    nkry = i
    niter = nkry+1

    ! 6) compute the coefficient of the Krylov expansion (back sub.)
    krylov: DO i = nkry, 1, -1
      ! coefficients are stored in c for convenience
      ci = rn2(i)
      DO k = nkry, i+1, -1
        ci = ci - h(i,k)*c(k)
      ENDDO
      c(i) = ci/h(i,i)
    ENDDO krylov

    ! 7) evaluate the Krylov expansion
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(jb, i) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = 1, no_of_blocks
      DO i = 1, nkry
        x(:,jb) = x(:,jb) + c(i)*v(:,jb,i)
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL 
    res(1) = ABS(rn2(niter))
    CALL trans%sync(x)
!    stop_timer(timer_gmres,2)

  END SUBROUTINE ocean_restart_gmres

END MODULE mo_ocean_solve_legacy_gmres
