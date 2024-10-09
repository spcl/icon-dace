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

! Estimate floating point error in RBF interpolation. This module
! helps to decide whether Cholesky applied to the RBF matrices will
! succeed in floating point arithmetic and gives a qualtiy measure
! for the chosen shape parameter.

MODULE mo_rbf_errana

  USE mo_kind,                 ONLY: wp
  USE mo_impl_constants,       ONLY: SUCCESS
  USE mo_math_constants,       ONLY: pi_180
  USE mo_exception,            ONLY: finish
  USE mo_parallel_config,      ONLY: nproma
  USE mo_model_domain,         ONLY: t_patch
  USE mo_communication,        ONLY: idx_1d
  USE mo_math_types,           ONLY: t_geographical_coordinates,                 &
    &                                t_cartesian_coordinates
  USE mo_math_utilities,       ONLY: gc2cc
  USE mo_math_utility_solvers, ONLY: solve_chol_v, solve_chol, choldec, choldec_v

  IMPLICIT NONE
  PUBLIC :: estimate_rbf_parameter

  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_rbf_errana')


CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Gaussian kernel for RBF interpolation.
  !!
  !! @f$\phi(r)=e^{-r^2}@f$
  !!
  FUNCTION gaussi (p_x, p_scale)  RESULT (p_rbf_val)
    REAL(wp) , INTENT(in) :: p_x             ! radial distance
    REAL(wp) , INTENT(in) :: p_scale         ! scale parameter
    REAL(wp)              :: p_rbf_val       ! RBF value

    p_rbf_val = p_x / p_scale
    p_rbf_val = -1._wp * p_rbf_val * p_rbf_val
    p_rbf_val = EXP(p_rbf_val)
  END FUNCTION gaussi


  ! -----------------------------------------------------------------
  !
  ! Estimate numerical error and interpolation error for radial
  ! basis functions (RBF).
  !
  ! INPUT:
  ! e     : RBF shape parameter
  ! xx,ff : data sites and function values
  !
  ! OUTPUT
  ! s     : stability value ([Demmel1989]):
  !         >0 : Cholesky decomposition is successful
  !
  SUBROUTINE rbf_error(e,start_idx,end_idx, kdim, max_kdim,jb, center, &
    &                  intp_data_iidx, intp_data_iblk, s, lflag)

    REAL(wp),              INTENT(IN)             :: e(:)                     !< shape parameter (1..nproma)
    INTEGER,               INTENT(IN)             :: start_idx,end_idx        !< start and end indices
    INTEGER,               INTENT(IN)             :: kdim(:)                  !< (nproma) matrix dimension
    INTEGER,               INTENT(IN)             :: max_kdim                 !< maximum matrix dimension
    INTEGER,               INTENT(IN)             :: jb                       !< block index
    TYPE (t_geographical_coordinates), INTENT(IN) :: center(:,:)              !< cell or edge center
    INTEGER,               INTENT(IN)             :: intp_data_iidx(:,:,:), &
      &                                              intp_data_iblk(:,:,:)    !< Indices of interpol. source points
    REAL(wp),              INTENT(INOUT)          :: s(:)                     !< stability(e) (1..nproma)
    LOGICAL,               INTENT(IN)             :: lflag(:)                 !< skip all indices with lflag==.TRUE.
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"rbf_error"
    REAL(wp), PARAMETER                         :: eps = 1.e-15_wp          !< machine epsilon
                                                
    REAL(wp), ALLOCATABLE                       :: mat_A(:,:,:), z_diag(:,:), e_n(:,:),   &
      &                                            sign_v(:,:), gamma(:), gamma_2(:),     &
      &                                            a(:), b(:), v(:,:), p(:,:), q(:,:)
    REAL(wp)                                       d1, d2
    INTEGER                                     :: N, j1, j2, p_imax(1), i, jc, iN,       &
      &                                            iidx, iblk, ierrstat
    TYPE (t_cartesian_coordinates), ALLOCATABLE :: cc(:,:)
    TYPE (t_cartesian_coordinates)              :: dfr
    TYPE (t_geographical_coordinates)           :: gc(nproma)

    N = max_kdim
    ALLOCATE(mat_A(nproma,N,N), z_diag(nproma,N), e_n(nproma,N),               &
      &      sign_v(nproma,N), gamma(nproma), gamma_2(nproma), cc(nproma, N),  &
      &      a(nproma), b(nproma), v(N,nproma), p(N,nproma), q(N,nproma),      &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed')

    s(start_idx:end_idx)         = 0._wp
    mat_A(start_idx:end_idx,:,:) = 0._wp
    z_diag(start_idx:end_idx,:)  = 0._wp

    ! case 1: cell-based interpolation stencil
    DO j1 = 1,N
      DO jc=start_idx,end_idx
        IF (lflag(jc)) CYCLE
        IF (j1 > kdim(jc)) CYCLE
        ! data site j1, cell jc,jb
        !
        ! get lon/lat coordinates on unit sphere
        ! with index/block of data site j1
        iidx   = intp_data_iidx(j1,jc,jb)
        iblk   = intp_data_iblk(j1,jc,jb)
        gc(jc)%lon = center(iidx, iblk)%lon
        gc(jc)%lat = center(iidx, iblk)%lat
        ! get Cartesian coordinates
        cc(jc,j1) = gc2cc(gc(jc))
      END DO
    END DO

    ! compute RBF interpolation matrix A, distance matrix D, and
    ! matrix R (intermediate result)
    !
    ! (symmetric matrix A in R^{N x N}, where we have N data sites)
    DO j1 = 1,N
      DO j2 = 1,N
        DO jc=start_idx,end_idx
          IF (lflag(jc)) CYCLE
          IF ((j1>kdim(jc)) .OR. (j2>kdim(jc))) CYCLE
          ! evaluate scalar RBF, fill entry in A matrix
          dfr%x(:) = cc(jc,j1)%x(:) - cc(jc,j2)%x(:)
          ! d1 = norm(xx(i,:) - xx(j,:), 2)
          d1       = SQRT( DOT_PRODUCT(dfr%x, dfr%x) )
          ! A = phi(D,e)
          d2       = gaussi(d1, 1._wp/e(jc))
          mat_A(jc,j1,j2) = d2
        END DO
      END DO
    END DO
    ! decompose linear system
#ifdef __SX__
!CDIR NOIEXPAND
    CALL choldec_v(start_idx,end_idx,kdim,N,mat_A,z_diag)
#else
    CALL choldec_v(start_idx,end_idx,kdim,  mat_A,z_diag)
#endif

    ! --------------------------------------------------------------------
    ! first, estimate ||A^{-1}||_1 with a single step of Hager's
    ! algorithm optimizing | A^-1 x |_1 over the convex set
    !  S = ( x: |x|_1 <= 1 ).
    ! --------------------------------------------------------------------

    ! choose boundary point of the convex set of feasible points
    e_n(start_idx:end_idx,:) = 1.0_wp
    v(:,start_idx:end_idx)   = 0.0_wp
    p(:,start_idx:end_idx)   = 0.0_wp

#ifdef __SX__
!CDIR NOIEXPAND
    CALL solve_chol_v(start_idx, end_idx, kdim, N, mat_A, z_diag, e_n, v)
#else
    CALL solve_chol_v(start_idx, end_idx, kdim,    mat_A, z_diag, e_n, v)
#endif
    sign_v(start_idx:end_idx,:) = 0._wp

    DO i=1,N
      DO jc=start_idx,end_idx
        IF (lflag(jc)) CYCLE
        IF (v(i,jc) > 0._wp) sign_v(jc,i) =  1.0_wp
        IF (v(i,jc) < 0._wp) sign_v(jc,i) = -1.0_wp
      END DO
    END DO
#ifdef __SX__
!CDIR NOIEXPAND
    CALL solve_chol_v(start_idx, end_idx, kdim, N, mat_A, z_diag, sign_v, p)
#else
    CALL solve_chol_v(start_idx, end_idx, kdim,    mat_A, z_diag, sign_v, p)
#endif
    ! move to a new point, maximizing ||A^(-1) p||_1
    e_n(start_idx:end_idx,:) = 0._wp
    DO jc=start_idx,end_idx
      IF (lflag(jc)) CYCLE
      p_imax = MAXLOC(ABS(p(:,jc)))
      e_n(jc, p_imax(1)) = 1.0_wp
    END DO

#ifdef __SX__
!CDIR NOIEXPAND
    CALL solve_chol_v(start_idx, end_idx, kdim, N, mat_A, z_diag, e_n, v)
#else
    CALL solve_chol_v(start_idx, end_idx, kdim,    mat_A, z_diag, e_n, v)
#endif
    gamma(start_idx:end_idx) = 0._wp
    DO i=1,N
      DO jc=start_idx,end_idx
        IF (lflag(jc)) CYCLE
        IF (i>kdim(jc)) CYCLE
        gamma(jc) = gamma(jc) + ABS(v(i,jc))
      END DO
    END DO

    a(start_idx:end_idx)     = -1.0_wp
    b(start_idx:end_idx)     =  1.0_wp
    e_n(start_idx:end_idx,:) =  0.0_wp
    DO i=1,N
      DO jc=start_idx,end_idx
        IF (lflag(jc)) CYCLE
        IF (i>kdim(jc)) CYCLE
        e_n(jc,i) = a(jc)*b(jc)
        a(jc)     = -1.0_wp * a(jc)
        b(jc)     = b(jc) + 1._wp/(kdim(jc)-1)
      END DO
    END DO
    ! solve q = A^(-1)*e_n
#ifdef __SX__
!CDIR NOIEXPAND
    CALL solve_chol_v(start_idx, end_idx, kdim, N, mat_A, z_diag, e_n, q)
#else
    CALL solve_chol_v(start_idx, end_idx, kdim,    mat_A, z_diag, e_n, q)
#endif
    ! counteract numerical cancellation [Higham1988]
    gamma_2(start_idx:end_idx) = 0._wp
    DO i=1,N
      DO jc=start_idx,end_idx
        IF (lflag(jc)) CYCLE
        IF (i>kdim(jc)) CYCLE
        gamma_2(jc) = gamma_2(jc) + ABS(q(i,jc))
      END DO
    END DO

    DO jc=start_idx,end_idx
      IF (lflag(jc)) CYCLE
      iN = kdim(jc)
      gamma(jc) = MAX(gamma(jc), 2./(3*iN) * gamma_2(jc))  ! estimate ||A^{-1}||_1

      ! --------------------------------------------------------------------
      ! based on the above estimate "gamma" for ||A^{-1}||_1 we apply
      ! Demmel's threshold, checking if we can still guarantee the Cholesky
      ! decomposition to be successful in floating point arithmetic:
      ! --------------------------------------------------------------------
      s(jc) = 1._wp/gamma(jc) - SQRT(REAL(iN,wp))*eps*REAL(iN,wp)*(REAL(iN,wp)+1._wp)/(1._wp-(REAL(iN,wp)-1._wp)*eps)
    END DO

    DEALLOCATE(mat_A, z_diag, e_n, sign_v, gamma, gamma_2,  &
      &        a, b, v, p, q, cc, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed')
  END SUBROUTINE rbf_error


  ! ----------------------------------------------------------------------------
  ! RBF shape parameter estimation
  ! ----------------------------------------------------------------------------
  !
  ! Based on an analysis of the round-off error in floating point
  ! arithmetic.
  !
  ! DESCRIPTION OF THE ALGORITHM
  !
  ! We assume Gaussian RBF shape functions of the form
  !   f(r) = exp(-(c*r).^2)
  ! with shape parameter c. Note that the shape parameter is often also
  ! defined as the reciprocal value of c. (For example, this is the
  ! case in ICON's main interpolation routines. Therefore this
  ! subroutine returns the inverse of the computed c!).
  !
  ! This routine estimates the smallest value for the shape parameter
  ! for which the Cholesky is likely to succeed in floating point
  ! arithmetic. For this estimate, we use the error analysis in
  !
  ! Demmel, J.: "On Floating Point Errors in Cholesky" (1989)
  !
  ! where we estimate the inverse matrix 1-norm by a single step of
  ! Hager's algorithm, see, e.g.,
  !
  ! Higham, N. J.: "FORTRAN codes for estimating the one-norm of a 
  !                 real or complex matrix, with applications to 
  !                 condition estimations."
  !                 ACM Trans. Math. Softw. , 1988, 14, 381-396
  !
  ! We assume that for sufficiently small values of c there exists a
  ! log-linear relation of the form
  !   log(t(c)) = q*log(c) + beta
  ! with coefficients q, beta.
  ! 
  ! Thus, we start with a sufficiently large value and compute
  ! Demmel's error threshold t(c) for a finite sequence of shape
  ! parameters
  !   c0, c_fak*c0, c_fak^2*c0, ..., c_fak^(n-1)*c0
  ! From this sequence we compute the coefficients q, beta by
  ! linear regression.
  ! 
  ! When the correlation coefficient is sufficiently close to 1
  ! (measured by a given threshold value tol_r), we stop this
  ! iteration, otherwise the first value in our sequence is discarded
  ! and a new one is calculated. However, we stop at a given lower
  ! threshold for c.
  ! 
  ! Having obtained an approximation for the function t(c), we can the
  ! estimate (extrapolate) the shape parameter which fulfils a given
  ! stability threshold tol_c1.
  !
  ! ----------------------------------------------------------------------------
  SUBROUTINE estimate_rbf_parameter(dst_nblks_c, dst_npromz_c, center,                    &
    &                              intp_data_iidx, intp_data_iblk, intp_data_nstencil,    &
    &                              max_nstencil, global_idx, rbf_shape_param)
    INTEGER,               INTENT(IN)           :: dst_nblks_c, dst_npromz_c  !< size of destination grid
    TYPE (t_geographical_coordinates), INTENT(IN) :: center(:,:)              !< cell or edge center
    INTEGER,               INTENT(IN)           :: intp_data_iidx(:,:,:), &
      &                                            intp_data_iblk(:,:,:), &   !< Indices of interpolation source points
      &                                            intp_data_nstencil(:,:)
    INTEGER,               INTENT(IN)           :: max_nstencil               !< max. stencil size
    INTEGER,               INTENT(IN)           :: global_idx(:)              !< for each lon-lat point: global idx 
    REAL(wp),              INTENT(INOUT)        :: rbf_shape_param

    ! constants, defining the behavior of the algorithm
    INTEGER,  PARAMETER :: min_tests        =         50
    INTEGER,  PARAMETER :: max_tests_nproma =       1000  ! max. no. of tested cell block lines (jb's)
    !
    REAL(wp), PARAMETER :: c0           =  1000.0_wp  ! initial VALUE for shape parameter search
    REAL(wp), PARAMETER :: c_fak        =     0.9_wp  ! factor for generating the sequence
    INTEGER,  PARAMETER :: n            =       5     ! length of the sequence
    REAL(wp), PARAMETER :: tol_r        =   1.e-3_wp  ! threshold for correlation coefficient
    REAL(wp), PARAMETER :: tol_c0       =      2._wp  ! lower threshold for shape parameter search
    REAL(wp), PARAMETER :: tol_c1       =  1.e-10_wp  ! given stability threshold
    REAL(wp), PARAMETER :: tol_resample =   1.e-7_wp  ! threshold for double-log regression: add new sample

    ! local variables
    REAL(wp) :: r(nproma), a1, a2, q(nproma), beta(nproma), c_tol(nproma), &
      &         sum_y, sum_z, denom
    REAL(wp) :: c_seq(nproma,n), t_seq(nproma,n), z(n), y(n), result_val(dst_nblks_c)
    INTEGER  :: i, start_idx, end_idx, jc, jb, itest_stride, kdim(nproma), max_tests, &
      &         min_stencil(1)
    LOGICAL  :: lflag(nproma)

    ! for regional grids, a process may have something to do...
    IF ((dst_nblks_c == 0) .OR. ((dst_nblks_c == 1) .AND. (dst_npromz_c == 0))) THEN
      rbf_shape_param = 99.
      RETURN
    END IF
    
    max_tests = MAX(min_tests,min_tests*max_tests_nproma/nproma)
    result_val(:) = 0._wp
    itest_stride = MAX(1, dst_nblks_c/max_tests)
!$OMP PARALLEL DO PRIVATE(jb, start_idx, end_idx, min_stencil, c_seq, t_seq, q, &
!$OMP                     beta, lflag, jc, kdim, r, z, sum_z, a2, y, sum_y, a1, denom)
    BLOCKS: DO jb = 1,dst_nblks_c
      start_idx = 1
      end_idx   = nproma
      if (jb == dst_nblks_c) end_idx = dst_npromz_c

      min_stencil = MINLOC(intp_data_nstencil(start_idx:end_idx,jb))

      c_seq(start_idx:end_idx,:) = 0._wp
      t_seq(start_idx:end_idx,:) = 0._wp
      c_seq(start_idx:end_idx,2) = c0

      ! Initialization of variables that are used in the non-flagged final loop
      q(:)    = tol_c1
      beta(:) = LOG(tol_c1)

      ! define control samples: skip all indices with lflag==.TRUE.
      ! 
      ! Note: Our choice here should be processor-independent!
      !
      lflag(:) = .TRUE.
      DO jc=start_idx,end_idx
        lflag(jc) = (MOD(global_idx(idx_1d(jc,jb)), itest_stride) /= 0)
      END DO
      ! add the cell with the smallest stencil in this block
      ! (note: this is not processor-independent!)
      lflag(min_stencil(1)) = .FALSE.

      kdim(:) = 0
      ! compute only for values in our control sample:
      kdim(start_idx:end_idx) = MERGE(0, intp_data_nstencil(start_idx:end_idx,jb), &
        &                             lflag(start_idx:end_idx)) 
      
      CALL rbf_error(c_seq(:,2),start_idx,end_idx,kdim,max_nstencil,jb, &
        &            center, intp_data_iidx, intp_data_iblk, t_seq(:,2), lflag)
      DO i=3,n
        c_seq(start_idx:end_idx,i) = MERGE(0._wp, c_fak * c_seq(start_idx:end_idx,i-1), &
          &                                lflag(start_idx:end_idx))
        CALL rbf_error(c_seq(:,i),start_idx,end_idx,kdim,max_nstencil,jb, &
          &            center, intp_data_iidx, intp_data_iblk, t_seq(:,i), lflag)
      END DO
      r(:) = 0._wp
      DO
        IF (ALL(lflag(start_idx:end_idx))) EXIT

        kdim(:) = 0
        ! compute only for values in our control sample:
        kdim(start_idx:end_idx) = MERGE(0, intp_data_nstencil(start_idx:end_idx,jb), &
          &                             lflag(start_idx:end_idx)) 

        ! add a new value to the sequence
        DO jc=start_idx,end_idx
          IF (lflag(jc)) CYCLE
          c_seq(jc,1:(n-1))  = c_seq(jc,2:n)
          c_seq(jc,n)        = c_fak * c_seq(jc,n-1)
          t_seq(jc,1:(n-1))  = t_seq(jc,2:n)
        END DO
        CALL rbf_error(c_seq(:,n),start_idx,end_idx,kdim,max_nstencil,jb, &
          &            center, intp_data_iidx, intp_data_iblk, t_seq(:,n), lflag)

        ! linear regression
        DO jc=start_idx,end_idx
          IF (lflag(jc)) CYCLE
          z        = LOG(c_seq(jc,:))
          sum_z    = SUM(z)
          a2       = (REAL(n,wp) * DOT_PRODUCT(z,z) - sum_z*sum_z)
          y        = LOG(t_seq(jc,:))
          sum_y    = SUM(y)
          a1       = (REAL(n,wp) * DOT_PRODUCT(y,z) - sum_y*sum_z)
          q(jc)    = a1/a2
          beta(jc) = (sum_y - q(jc)*sum_z)/n
          ! In exact arithmetics the following check would not be
          ! necessary due to the Cauchy Schwarz inequality:
          IF ((REAL(n,wp)*DOT_PRODUCT(y,y) - sum_y*sum_y) < 0._wp) THEN
            CYCLE
          END IF
          denom    = SQRT(a2)*SQRT(REAL(n,wp)*DOT_PRODUCT(y,y) - sum_y*sum_y)
          IF (denom < tol_resample) THEN
            ! now, there seems to be a problem: the function samples
            ! do not reveal a double-log form.
            !
            ! skip the current test and add another (smaller) test
            ! value to our set of samples:
            CYCLE
          END IF
          r(jc)    = a1/denom  ! correlation coefficient
          IF (c_seq(jc,n)  <= tol_c0) lflag(jc) = .TRUE.
          IF (ABS(1-r(jc)) <= tol_r)  lflag(jc) = .TRUE.
        END DO
      END DO
      ! compute the shape parameter which fulfils a given stability
      ! threshold:
      c_tol(start_idx:end_idx) = 0._wp
      DO jc=start_idx,end_idx
        c_tol(jc) = MAX(tol_c0, EXP((LOG(tol_c1) - beta(jc))/q(jc)))
      END DO
      result_val(jb) = MAXVAL(c_tol(start_idx:end_idx))
    END DO BLOCKS
!$OMP END PARALLEL DO
    rbf_shape_param = 1._wp/MAXVAL(result_val)
  END SUBROUTINE estimate_rbf_parameter
  
END MODULE mo_rbf_errana
