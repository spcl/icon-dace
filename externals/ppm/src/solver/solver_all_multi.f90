!>
!! @file solver_all_multi.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!! @version 1.0
!
! Keywords: scales ppm solver
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
#define filename 'solver_all_multi.f90'
!
! Dispatcher solve function which calls the appropiate solver/preconditioner
! according to the config settings
FUNCTION SOLVE(A, b, x, ext_x, exchange, tol_opt, maxiter_opt) RESULT(kiter)
  USE solver_internal
  USE ppm_extents, ONLY: extent
  USE ppm_base, ONLY: abort_ppm

  REAL(PREC), INTENT(IN) :: b(:,:)                ! right-hand-side b
  REAL(PREC), INTENT(INOUT) :: x(:,:)             ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(PREC), OPTIONAL, INTENT(IN) :: tol_opt     ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
  INTEGER :: kiter, maxiter
  REAL(PREC) :: tol, lambda_min, lambda_max

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
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

  kiter = 0
  ! Check and set optional arguments
  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt
  ELSE
    tol = REAL(config%tol, PREC)
  ENDIF

  ! Make sure stencil is set up
  IF (.NOT. stencil_defined(0.0_precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         filename, __LINE__)
  ENDIF

  ! Make sure preconditioner is initialized
  IF (.NOT. precond_prepared(config%preconditioner, 0.0_precs)) THEN
    SELECT CASE (config%preconditioner)
    CASE (JACOBI_PRECOND)
      CALL prep_jacobi(0.0_precs)
    CASE (ILU0_PRECOND)
      CALL prep_ilu0(0.0_precs)
    CASE (SSOR_PRECOND)
      CALL abort_ppm("No SSOR parameter provided!", filename, __LINE__)
    CASE (ICC_PRECOND)
      CALL prep_icc(config%icc_param, 0.0_precs)
    CASE (MICC_PRECOND)
      CALL prep_micc(config%icc_param, 0.0_precs)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "]  does not exist!", filename, __LINE__)
    END SELECT
  ENDIF

  ! Call right preconditioner/solver combination
  SELECT CASE (config%solver)
  CASE (CG_SOLVER)
    SELECT CASE (config%preconditioner)
    CASE (NONE_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, IDENTITY, &
           exchange, tol, maxiter)
    CASE (JACOBI_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, JACOBI, &
           exchange, tol, maxiter)
    CASE (ILU0_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, ILU0, &
           exchange, tol, maxiter)
    CASE (SSOR_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, SSOR, &
           exchange, tol, maxiter)
    CASE (ICC_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, ICC, &
           exchange, tol, maxiter)
    CASE (MICC_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, ICC, &
           exchange, tol, maxiter)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "] does not exist!", filename, __LINE__)
    END SELECT
  CASE (CHEBYSHEV_SOLVER)
    lambda_min = REAL(config%lambda_min, PREC)
    lambda_max = REAL(config%lambda_max, PREC)
    IF (ieee_is_nan(lambda_min) .OR. ieee_is_nan(lambda_max)) THEN
      CALL abort_ppm("Lambda_min and lambda_max not provided!", filename, &
           __LINE__)
    ENDIF
    SELECT CASE (config%preconditioner)
    CASE (NONE_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, IDENTITY, &
           exchange, tol, maxiter)
    CASE (JACOBI_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, JACOBI, &
           exchange, tol, maxiter)
    CASE (ILU0_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, ILU0, &
           exchange, tol, maxiter)
    CASE (SSOR_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, SSOR, &
           exchange, tol, maxiter)
    CASE (ICC_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, ICC, &
           exchange, tol, maxiter)
    CASE (MICC_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, ICC, &
           exchange, tol, maxiter)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "] does not exist!", filename, __LINE__)
    END SELECT
  CASE DEFAULT
    CALL abort_ppm("Selected solver [" // int2str(config%solver) &
         // "]  does not exist!", filename, __LINE__)
  END SELECT

END FUNCTION SOLVE

#undef filename
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
