!>
!! @file spectral_methods.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!! @author Florian Wilhelm <Florian.Wilhelm@kit.edu>
!
! Keywords: scales ppm solver eigenvalue sor
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
!
!> provides functions to determine eigenvalues

MODULE spectral_methods
  USE solver_internal
  USE linear_algebra
  USE ppm_std_type_kinds, ONLY: sp, dp
  USE ppm_extents, ONLY: extent, extent_start, extent_end
  USE ppm_base, ONLY: abort_ppm
  USE preconditioners, ONLY: precond_prepared, jacobi_precond, &
       prep_jacobi,jacobi
  USE solver_config, ONLY: stencil_defined

  IMPLICIT NONE

  PRIVATE

  !-------------------------------------------------------------------
  ! Module Function & Procedure Interface
  !-------------------------------------------------------------------

  !> @brief power method to determine largest absolute eigenvalue and vector
  !! @details variants in single and double precision
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in,out] x start vector, holds eigenvector in the end
  !! @param[in] ext_x extent of x
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return eigenvalue
  INTERFACE power_method
    MODULE PROCEDURE power_method_sp
    MODULE PROCEDURE power_method_dp
  END INTERFACE

  !> @brief calculates largest absolute eigenvalue
  !! @details variants in single and double precision, same as power method
  !! but leaves x unchanged
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] x start vector
  !! @param[in] ext_x extent of x
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return largest absolute eigenvalue
  INTERFACE calc_lambda_max
    MODULE PROCEDURE calc_lambda_max_sp
    MODULE PROCEDURE calc_lambda_max_dp
  END INTERFACE

  !> @brief calculates smallest absolute eigenvalue
  !! @details variants in single and double precision
  !! @param[in] A_shift shifted @link solver_internal::linop linear operator
  !! @endlink like matrix multiplication
  !! @param[in] x start vector
  !! @param[in] ext_x extent of x
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] lambda_max largest absolute eigenvalue
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return smallest absolute eigenvalue
  INTERFACE calc_lambda_min
    MODULE PROCEDURE calc_lambda_min_sp
    MODULE PROCEDURE calc_lambda_min_dp
  END INTERFACE

  !> @brief calculates smallest and largest absolute eigenvalue
  !! @details variants in single and double precision
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] A_shift shifted @link solver_internal::linop linear operator
  !! @endlink like matrix multiplication
  !! @param[in] x start vector
  !! @param[in] ext_x extent of x
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return array with (lambda_min, lambda_max)
  INTERFACE calc_lambda_min_max
    MODULE PROCEDURE calc_lambda_min_max_sp
    MODULE PROCEDURE calc_lambda_min_max_dp
  END INTERFACE

  !> @brief calculates optimal SOR parameter
  !! @details variants in single and double precision
  !! @param[in] rho spectral radius of Jacobi iteration matrix
  !! @return SOR parameter
  INTERFACE calc_optimal_sor_param
    MODULE PROCEDURE calc_optimal_sor_param_sp
    MODULE PROCEDURE calc_optimal_sor_param_dp
  END INTERFACE

  !> @brief Jacobi iteration matrix
  !! @param[in] x 2d input array
  !! @param[out] y 2d output array
  INTERFACE jacobi_iter
    MODULE PROCEDURE jacobi_iter_sp
    MODULE PROCEDURE jacobi_iter_dp
  END INTERFACE

  PUBLIC:: power_method, calc_lambda_max, calc_lambda_min &
       , calc_lambda_min_max, calc_optimal_sor_param, jacobi_iter &
       , jacobi_iter_sp, jacobi_iter_dp

CONTAINS

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
!
! Power method to calculate dominant eigenvalue and corresponding eigenvector
FUNCTION power_method_sp(A, x, ext_x, exchange, tol_opt, maxiter_opt) &
     RESULT(lambda_max)
  REAL(sp), INTENT(INOUT) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(sp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt    ! maximum number of iterations
  REAL(sp), ALLOCATABLE :: y(:,:)
  REAL(sp) :: lambda_max, lambda_old, lambda_new, tol, norm
  INTEGER :: kiter, maxiter, ib, il, jb, jl, ie, je

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(IN) :: field(:,:)
      REAL(sp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, sp)
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

END FUNCTION power_method_sp

! Same as powermethod but leaves x unchanged, calculates dominant eigenvalue
FUNCTION calc_lambda_max_sp(A, x, ext_x, exchange, tol_opt, maxiter_opt) &
     RESULT(lambda_max)
  REAL(sp), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(sp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(sp), ALLOCATABLE :: y(:,:)
  REAL(sp) :: lambda_max, tol
  INTEGER :: maxiter, ie, je

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(IN) :: field(:,:)
      REAL(sp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, sp)
  ENDIF

  y = x

  lambda_max = power_method(A, y, ext_x, exchange, tol, maxiter)

  DEALLOCATE(y)

END FUNCTION calc_lambda_max_sp

! Calculates least-dominant eigenvalue with the help of the dominant eigenvalue
FUNCTION calc_lambda_min_sp(A_shifted, x, ext_x, exchange, lambda_max, tol_opt, &
     maxiter_opt) RESULT(lambda_min)
  REAL(sp), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(sp), INTENT(IN) :: lambda_max
  REAL(sp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(sp) :: lambda_min, tol
  INTEGER :: maxiter

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  INTERFACE
    ! Matrix-Vector multiplication of the linear system shifted by a value
    SUBROUTINE A_shifted(field, res_field)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(IN) :: field(:,:)
      REAL(sp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A_shifted
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, sp)
  ENDIF

  ! Shift by stencil by lambda_max
  stencil_sp%shift = lambda_max

  ! Calculate dominant eigenvalue with altered stencil
  lambda_min = calc_lambda_max(A_shifted, x, ext_x, exchange, tol, maxiter)
  lambda_min = lambda_max - lambda_min

  ! Reset stencil shift
  stencil_sp%shift = 0.0_precs

END FUNCTION calc_lambda_min_sp

! Calculates dominant and least dominant eigenvalues
FUNCTION calc_lambda_min_max_sp(A, A_shifted, x, ext_x, exchange, tol_opt, &
     maxiter_opt) RESULT(lambdas)
  REAL(sp), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(sp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(sp) :: tol
  REAL(sp), DIMENSION(2) :: lambdas
  INTEGER :: maxiter

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(IN) :: field(:,:)
      REAL(sp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Matrix-Vector multiplication of the linear system shifted by a value
    SUBROUTINE A_shifted(field, res_field)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(IN) :: field(:,:)
      REAL(sp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A_shifted
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: sp
      REAL(sp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, sp)
  ENDIF

  lambdas(2) = calc_lambda_max(A, x, ext_x, exchange, tol, maxiter)
  lambdas(1) = calc_lambda_min(A_shifted, x, ext_x, exchange, lambdas(2), tol, &
       maxiter)

END FUNCTION calc_lambda_min_max_sp

! Calculate optimal SOR parameter which can also be used as SSOR parameter
! Attention: Keep in mind that the spectral radius rho is the one
!            of the jacobi iteration matrix B = D^-1*(-L -U) [A = D + L + U]
!            See [SAAD, p.105] for references
FUNCTION calc_optimal_sor_param_sp(rho) RESULT(sor_param)
  REAL(sp), INTENT(IN) :: rho
  REAL(sp) :: sor_param

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  IF ( rho < 0.0_precs .OR. rho >= 1.0_precs ) THEN
    CALL abort_ppm("Provided spectral radius rho does not satisfy 0<=rho<1.", &
         'spectral_methods_multi.f90', 298)
  ENDIF

  sor_param = 2.0_precs/(1.0_precs + SQRT(1.0_precs-rho**2.0_precs))

END FUNCTION calc_optimal_sor_param_sp

! The jacobi iteration matrix B as needed to calculate optimal SOR parameter
! B = D^-1*(-L -U) [A = D + L + U] with D=diagonal, L=lower, U=upper
SUBROUTINE jacobi_iter_sp(field, res_field)

  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: ib, il, jb, jl, i, j
  REAL(sp), DIMENSION(:,:), POINTER :: uf, vf

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  IF (.NOT. stencil_defined(0._precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         'spectral_methods_multi.f90', 319)
  ENDIF

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Define some aliases
  uf => stencil_sp%zonal
  vf => stencil_sp%meridional

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

END SUBROUTINE jacobi_iter_sp

!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
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
!
! Power method to calculate dominant eigenvalue and corresponding eigenvector
FUNCTION power_method_dp(A, x, ext_x, exchange, tol_opt, maxiter_opt) &
     RESULT(lambda_max)
  REAL(dp), INTENT(INOUT) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(dp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt    ! maximum number of iterations
  REAL(dp), ALLOCATABLE :: y(:,:)
  REAL(dp) :: lambda_max, lambda_old, lambda_new, tol, norm
  INTEGER :: kiter, maxiter, ib, il, jb, jl, ie, je

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(IN) :: field(:,:)
      REAL(dp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, dp)
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

END FUNCTION power_method_dp

! Same as powermethod but leaves x unchanged, calculates dominant eigenvalue
FUNCTION calc_lambda_max_dp(A, x, ext_x, exchange, tol_opt, maxiter_opt) &
     RESULT(lambda_max)
  REAL(dp), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(dp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(dp), ALLOCATABLE :: y(:,:)
  REAL(dp) :: lambda_max, tol
  INTEGER :: maxiter, ie, je

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(IN) :: field(:,:)
      REAL(dp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, dp)
  ENDIF

  y = x

  lambda_max = power_method(A, y, ext_x, exchange, tol, maxiter)

  DEALLOCATE(y)

END FUNCTION calc_lambda_max_dp

! Calculates least-dominant eigenvalue with the help of the dominant eigenvalue
FUNCTION calc_lambda_min_dp(A_shifted, x, ext_x, exchange, lambda_max, tol_opt, &
     maxiter_opt) RESULT(lambda_min)
  REAL(dp), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(dp), INTENT(IN) :: lambda_max
  REAL(dp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(dp) :: lambda_min, tol
  INTEGER :: maxiter

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  INTERFACE
    ! Matrix-Vector multiplication of the linear system shifted by a value
    SUBROUTINE A_shifted(field, res_field)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(IN) :: field(:,:)
      REAL(dp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A_shifted
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, dp)
  ENDIF

  ! Shift by stencil by lambda_max
  stencil_dp%shift = lambda_max

  ! Calculate dominant eigenvalue with altered stencil
  lambda_min = calc_lambda_max(A_shifted, x, ext_x, exchange, tol, maxiter)
  lambda_min = lambda_max - lambda_min

  ! Reset stencil shift
  stencil_dp%shift = 0.0_precs

END FUNCTION calc_lambda_min_dp

! Calculates dominant and least dominant eigenvalues
FUNCTION calc_lambda_min_max_dp(A, A_shifted, x, ext_x, exchange, tol_opt, &
     maxiter_opt) RESULT(lambdas)
  REAL(dp), INTENT(IN) :: x(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)         ! extent of x
  REAL(dp), INTENT(IN), OPTIONAL :: tol_opt
  INTEGER, INTENT(IN), OPTIONAL :: maxiter_opt ! maximum number of iterations
  REAL(dp) :: tol
  REAL(dp), DIMENSION(2) :: lambdas
  INTEGER :: maxiter

  INTERFACE
    ! Matrix-Vector multiplication of the linear system
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(IN) :: field(:,:)
      REAL(dp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Matrix-Vector multiplication of the linear system shifted by a value
    SUBROUTINE A_shifted(field, res_field)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(IN) :: field(:,:)
      REAL(dp), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A_shifted
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: dp
      REAL(dp), INTENT(INOUT) :: a0(:,:)
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
    tol = REAL(config%lambda_tol, dp)
  ENDIF

  lambdas(2) = calc_lambda_max(A, x, ext_x, exchange, tol, maxiter)
  lambdas(1) = calc_lambda_min(A_shifted, x, ext_x, exchange, lambdas(2), tol, &
       maxiter)

END FUNCTION calc_lambda_min_max_dp

! Calculate optimal SOR parameter which can also be used as SSOR parameter
! Attention: Keep in mind that the spectral radius rho is the one
!            of the jacobi iteration matrix B = D^-1*(-L -U) [A = D + L + U]
!            See [SAAD, p.105] for references
FUNCTION calc_optimal_sor_param_dp(rho) RESULT(sor_param)
  REAL(dp), INTENT(IN) :: rho
  REAL(dp) :: sor_param

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  IF ( rho < 0.0_precs .OR. rho >= 1.0_precs ) THEN
    CALL abort_ppm("Provided spectral radius rho does not satisfy 0<=rho<1.", &
         'spectral_methods_multi.f90', 298)
  ENDIF

  sor_param = 2.0_precs/(1.0_precs + SQRT(1.0_precs-rho**2.0_precs))

END FUNCTION calc_optimal_sor_param_dp

! The jacobi iteration matrix B as needed to calculate optimal SOR parameter
! B = D^-1*(-L -U) [A = D + L + U] with D=diagonal, L=lower, U=upper
SUBROUTINE jacobi_iter_dp(field, res_field)

  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: ib, il, jb, jl, i, j
  REAL(dp), DIMENSION(:,:), POINTER :: uf, vf

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  IF (.NOT. stencil_defined(0._precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         'spectral_methods_multi.f90', 319)
  ENDIF

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Define some aliases
  uf => stencil_dp%zonal
  vf => stencil_dp%meridional

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

END SUBROUTINE jacobi_iter_dp

!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:

END MODULE spectral_methods
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
