!>
!! @file spectral_methods_multi.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!
! Keywords: scales ppm solver eigenvalue sor
! Author: Florian Wilhelm <Florian.Wilhelm@kit.edu>
! Maintainer: Florian Wilhelm <Florian.Wilhelm@kit.edu>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
#define filename 'spectral_methods_multi.f90'
!
! Power method to calculate dominant eigenvalue and corresponding eigenvector
FUNCTION POWER_METHOD(A, x, ext_x, exchange, tol_opt, maxiter_opt) &
     RESULT(lambda_max)
  REAL(PREC), INTENT(INOUT) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(PREC), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt    ! maximum number of iterations
  REAL(PREC), ALLOCATABLE :: y(:,:)
  REAL(PREC) :: lambda_max, lambda_old, lambda_new, tol, norm
  INTEGER :: kiter, maxiter, ib, il, jb, jl, ie, je

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  ! Check and set optional arguments
  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%lambda_maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt
  ELSE
    tol = REAL(config%lambda_tol, PREC)
  ENDIF

  ie = SIZE(x,1)
  je = SIZE(x,2)
  ib = extent_start(ext_x(1))
  jb = extent_start(ext_x(2))
  il = extent_end(ext_x(1))
  jl = extent_end(ext_x(2))

  ALLOCATE(y(ie,je))

  ! Make sure y is ordered (not NaN nor +-Inf) outside of ib:il,jb:jl
  y(1:ib-1,:) = 0.0_precs
  y(il+1:ie,:) = 0.0_precs
  y(ib:il,1:jb-1) = 0.0_precs
  y(ib:il,jl+1:je) = 0.0_precs

  kiter = 0
  lambda_old = 1.0_precs
  lambda_new = 10.0_precs
  lambda_max = 0.0_precs

  norm = arr_norm_2(x(ib:il,jb:jl))
  IF (norm > 0.0_precs) THEN
    x(ib:il,jb:jl) = x(ib:il,jb:jl)/norm

    IF (config%do_exchange) CALL exchange(x, 'power_method 1')

    DO WHILE ( ABS(lambda_new - lambda_old) > tol .AND. kiter < maxiter )
      CALL A(x, y)  ! y = A*x
      norm = arr_norm_2(y(ib:il,jb:jl))
      IF (norm > 0.0_precs) THEN
        x(ib:il,jb:jl) = y(ib:il,jb:jl)/norm
        lambda_old = lambda_new
        lambda_new = arr_dotproduct(x(ib:il,jb:jl),y(ib:il,jb:jl))
        IF (config%do_exchange) CALL exchange(x, 'power_method 2')
        kiter = kiter + 1
        !WRITE(*,*) "kiter: ", kiter, "lambda_new: ", lambda_new
      ELSE
        RETURN
      END IF
    ENDDO

    lambda_max = lambda_new

  END IF

END FUNCTION POWER_METHOD

! Same as powermethod but leaves x unchanged, calculates dominant eigenvalue
FUNCTION CALC_LAMBDA_MAX(A, x, ext_x, exchange, tol_opt, maxiter_opt) &
     RESULT(lambda_max)
  REAL(PREC), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(PREC), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(PREC), ALLOCATABLE :: y(:,:)
  REAL(PREC) :: lambda_max, tol
  INTEGER :: maxiter, ie, je

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  ie = SIZE(x,1)
  je = SIZE(x,2)

  ALLOCATE(y(ie,je))

  ! Check and set optional arguments
  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%lambda_maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt
  ELSE
    tol = REAL(config%lambda_tol, PREC)
  ENDIF

  y = x

  lambda_max = power_method(A, y, ext_x, exchange, tol, maxiter)

  DEALLOCATE(y)

END FUNCTION CALC_LAMBDA_MAX

! Calculates least-dominant eigenvalue with the help of the dominant eigenvalue
FUNCTION CALC_LAMBDA_MIN(A_shifted, x, ext_x, exchange, lambda_max, tol_opt, &
     maxiter_opt) RESULT(lambda_min)
  REAL(PREC), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(PREC), INTENT(IN) :: lambda_max
  REAL(PREC), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(PREC) :: lambda_min, tol
  INTEGER :: maxiter

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  INTERFACE
    ! Matrix-Vector multiplication of the linear system shifted by a value
    SUBROUTINE A_shifted(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A_shifted
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  ! Check and set optional arguments
  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%lambda_maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt
  ELSE
    tol = REAL(config%lambda_tol, PREC)
  ENDIF

  ! Shift by stencil by lambda_max
  STENCIL%shift = lambda_max

  ! Calculate dominant eigenvalue with altered stencil
  lambda_min = calc_lambda_max(A_shifted, x, ext_x, exchange, tol, maxiter)
  lambda_min = lambda_max - lambda_min

  ! Reset stencil shift
  STENCIL%shift = 0.0_precs

END FUNCTION CALC_LAMBDA_MIN

! Calculates dominant and least dominant eigenvalues
FUNCTION CALC_LAMBDA_MIN_MAX(A, A_shifted, x, ext_x, exchange, tol_opt, &
     maxiter_opt) RESULT(lambdas)
  REAL(PREC), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(PREC), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(PREC) :: tol
  REAL(PREC), DIMENSION(2) :: lambdas
  INTEGER :: maxiter

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Matrix-Vector multiplication of the linear system shifted by a value
    SUBROUTINE A_shifted(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A_shifted
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  ! Check and set optional arguments
  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%lambda_maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt
  ELSE
    tol = REAL(config%lambda_tol, PREC)
  ENDIF

  lambdas(2) = calc_lambda_max(A, x, ext_x, exchange, tol, maxiter)
  lambdas(1) = calc_lambda_min(A_shifted, x, ext_x, exchange, lambdas(2), tol, &
       maxiter)

END FUNCTION CALC_LAMBDA_MIN_MAX

! Calculate optimal SOR parameter which can also be used as SSOR parameter
! Attention: Keep in mind that the spectral radius rho is the one
!            of the jacobi iteration matrix B = D^-1*(-L -U) [A = D + L + U]
!            See [SAAD, p.105] for references
FUNCTION CALC_OPTIMAL_SOR_PARAM(rho) RESULT(sor_param)
  REAL(PREC), INTENT(IN) :: rho
  REAL(PREC) :: sor_param

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  IF ( rho < 0.0_precs .OR. rho >= 1.0_precs ) THEN
    CALL abort_ppm("Provided spectral radius rho does not satisfy 0<=rho<1.", &
         filename, __LINE__)
  ENDIF

  sor_param = 2.0_precs/(1.0_precs + SQRT(1.0_precs-rho**2.0_precs))

END FUNCTION CALC_OPTIMAL_SOR_PARAM

! The jacobi iteration matrix B as needed to calculate optimal SOR parameter
! B = D^-1*(-L -U) [A = D + L + U] with D=diagonal, L=lower, U=upper
SUBROUTINE JACOBI_ITER(field, res_field)

  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)
  INTEGER :: ib, il, jb, jl, i, j
  REAL(PREC), DIMENSION(:,:), POINTER :: uf, vf

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  IF (.NOT. stencil_defined(0._precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         filename, __LINE__)
  ENDIF

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Define some aliases
  uf => STENCIL%zonal
  vf => STENCIL%meridional

  ! Apply matrix stencil with (-L -U)
  DO j=jb,jl
    DO i=ib,il
      res_field(i,j) = - uf(i,j)*field(i+1,j) - uf(i-1,j)*field(i-1,j) &
           - vf(i,j)*field(i,j+1) - vf(i,j-1)*field(i,j-1)
    ENDDO
  ENDDO

  IF (.NOT. precond_prepared(JACOBI_PRECOND, 0._precs)) THEN
    CALL prep_jacobi(0._precs)
  ENDIF

  ! Multiply with D^-1
  CALL jacobi(res_field)

END SUBROUTINE JACOBI_ITER

#undef filename
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
