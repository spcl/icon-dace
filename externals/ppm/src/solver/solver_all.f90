!>
!! @file solver_all.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
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
!
!> generic solve function and multi-precision iterative refinement

MODULE solver_all
  USE solver_config
  USE linear_algebra
  USE solvers
  USE preconditioners
  USE spectral_methods
  USE ieee_arithmetic, ONLY: ieee_is_nan
  USE ppm_std_type_kinds, ONLY: wp, sp, dp

  IMPLICIT NONE

  !> @brief general solve function to solve Ax = b in single/double precision
  !! @details dispatcher solve function which calls the appropiate
  !! solver/preconditioner according to the config settings
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] b right hand side
  !! @param[in,out] x start vector, holds solution in the end
  !! @param[in] ext_x extent of x
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return number of iterations
  INTERFACE solve
    MODULE PROCEDURE solve_sp
    MODULE PROCEDURE solve_dp
  END INTERFACE

CONTAINS

  !> @brief multi-precision iterative refinement
  !! @param[in] A_dp @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication in double precision
  !! @param[in] A_sp @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication in single precision
  !! @param[in] b right hand side in double precision
  !! @param[in,out] x start vector in double precision, hold solution in the end
  !! @param[in] ext_x extent of x
  !! @param[in] exchange_sp single precision @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return number of iterations
  FUNCTION iter_refinement(A_dp, A_sp, b, x, ext_x, exchange_sp, tol_opt, &
       maxiter_opt) RESULT(totaliter)
    USE solver_internal
    USE ppm_extents, ONLY: extent_start, extent_end, extent



    REAL(dp), INTENT(IN) :: b(:,:)                  ! right-hand-side b
    REAL(dp), INTENT(INOUT) :: x(:,:)               ! startvalue and result
    TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
    REAL(dp), OPTIONAL, INTENT(IN) :: tol_opt       ! tolerance for residual
    INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
    INTEGER :: jb, ib, il, jl, ie, je, totaliter, kiter, kkiter, maxiter



    REAL(dp) :: delta, delta0, tol, relres
    REAL(sp) :: slvprec
    REAL(dp), ALLOCATABLE :: r_dp(:,:), c_dp(:,:), Ax(:,:)
    REAL(sp), ALLOCATABLE :: r_sp(:,:), c_sp(:,:)

    INTERFACE
      ! Douple Prec. Matrix-Vector multiplication of the linear system to solve
      SUBROUTINE A_dp(field, res_field)
        USE ppm_std_type_kinds, ONLY: dp
        REAL(dp), INTENT(IN) :: field(:,:)
        REAL(dp), INTENT(OUT) :: res_field(:,:)
      END SUBROUTINE A_dp
      ! Single Prec. Matrix-Vector multiplication of the linear system to solve
      SUBROUTINE A_sp(field, res_field)
        USE ppm_std_type_kinds, ONLY: sp
        REAL(sp), INTENT(IN) :: field(:,:)
        REAL(sp), INTENT(OUT) :: res_field(:,:)
      END SUBROUTINE A_sp
      ! Function to exchange boundaries if necessary
      SUBROUTINE exchange_sp(a0, text)
        USE ppm_std_type_kinds, ONLY: sp
        REAL(sp), INTENT(INOUT) :: a0(:,:)
        CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
      END SUBROUTINE exchange_sp
    END INTERFACE

    ! Check and set optional arguments
    IF ( PRESENT(maxiter_opt) ) THEN
      maxiter = maxiter_opt
    ELSE
      maxiter = config%iterref_maxiter
    ENDIF
    IF ( PRESENT(tol_opt) ) THEN
      tol = tol_opt**2 ! **2, because we compare later the square of a norm
    ELSE
      tol = config%iterref_tol**2
    ENDIF

    ! Define some variables for the ranges of the fields
    ie = SIZE(x,1)
    je = SIZE(x,2)
    ib = extent_start(ext_x(1))
    jb = extent_start(ext_x(2))
    il = extent_end(ext_x(1))
    jl = extent_end(ext_x(2))




    ALLOCATE(r_dp(ie,je), r_sp(ie,je), c_dp(ie,je), c_sp(ie,je), Ax(ie,je))

    ! Initialise residual and direction
    CALL A_dp(x, Ax)

    r_dp(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)



    CALL clear_halos(r_sp, ext_x)

    ! Calculate scalar product of residual and direction
    delta0 = 1._dp/arr_dotproduct(r_dp(ib:il,jb:jl))

    relres = 1._dp
    kiter = 0
    totaliter = 0

    DO WHILE (kiter < maxiter .AND. relres > tol)
      ! tyepcast to low precision
      r_sp(ib:il,jb:jl) = REAL(r_dp(ib:il,jb:jl), sp)
      ! call bounds exchange
      CALL exchange_sp(r_sp, 'iter_refinement 1')
      ! Start value for single precision solver
      c_sp = 0._sp
      ! solve in low precision
      slvprec = REAL(MAX(config%tol,tol/relres), sp)
      kkiter = solve(A_sp, r_sp, c_sp, ext_x, exchange_sp, slvprec)
      ! typecast to high precision on the whole fields
      c_dp = REAL(c_sp, dp)
      ! perform the error correction
      x = x + c_dp
      ! calculate new residual
      CALL A_dp(x, Ax)

      r_dp(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)




      delta = arr_dotproduct(r_dp(ib:il,jb:jl))
      relres = delta*delta0
      kiter = kiter + 1
      totaliter = totaliter + kkiter

    ENDDO

    DEALLOCATE(r_dp, r_sp, c_dp, c_sp, Ax)

  END FUNCTION iter_refinement

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
!
! Dispatcher solve function which calls the appropiate solver/preconditioner
! according to the config settings
FUNCTION solve_sp(A, b, x, ext_x, exchange, tol_opt, maxiter_opt) RESULT(kiter)
  USE solver_internal
  USE ppm_extents, ONLY: extent
  USE ppm_base, ONLY: abort_ppm

  REAL(sp), INTENT(IN) :: b(:,:)                ! right-hand-side b
  REAL(sp), INTENT(INOUT) :: x(:,:)             ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(sp), OPTIONAL, INTENT(IN) :: tol_opt     ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
  INTEGER :: kiter, maxiter
  REAL(sp) :: tol, lambda_min, lambda_max

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
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
    tol = REAL(config%tol, sp)
  ENDIF

  ! Make sure stencil is set up
  IF (.NOT. stencil_defined(0.0_precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         'solver_all_multi.f90', 91)
  ENDIF

  ! Make sure preconditioner is initialized
  IF (.NOT. precond_prepared(config%preconditioner, 0.0_precs)) THEN
    SELECT CASE (config%preconditioner)
    CASE (JACOBI_PRECOND)
      CALL prep_jacobi(0.0_precs)
    CASE (ILU0_PRECOND)
      CALL prep_ilu0(0.0_precs)
    CASE (SSOR_PRECOND)
      CALL abort_ppm("No SSOR parameter provided!", 'solver_all_multi.f90', 102)
    CASE (ICC_PRECOND)
      CALL prep_icc(config%icc_param, 0.0_precs)
    CASE (MICC_PRECOND)
      CALL prep_micc(config%icc_param, 0.0_precs)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "]  does not exist!", 'solver_all_multi.f90', 110)
    END SELECT
  ENDIF

  ! Call right preconditioner/solver combination
  SELECT CASE (config%solver)
  CASE (CG_SOLVER)
    SELECT CASE (config%preconditioner)
    CASE (NONE_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, identity_sp, &
           exchange, tol, maxiter)
    CASE (JACOBI_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, jacobi_sp, &
           exchange, tol, maxiter)
    CASE (ILU0_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, ilu0_sp, &
           exchange, tol, maxiter)
    CASE (SSOR_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, SSOR_sp, &
           exchange, tol, maxiter)
    CASE (ICC_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, icc_sp, &
           exchange, tol, maxiter)
    CASE (MICC_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, icc_sp, &
           exchange, tol, maxiter)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "] does not exist!", 'solver_all_multi.f90', 139)
    END SELECT
  CASE (CHEBYSHEV_SOLVER)
    lambda_min = REAL(config%lambda_min, sp)
    lambda_max = REAL(config%lambda_max, sp)
    IF (ieee_is_nan(lambda_min) .OR. ieee_is_nan(lambda_max)) THEN
      CALL abort_ppm("Lambda_min and lambda_max not provided!", 'solver_all_multi.f90', &
           146)
    ENDIF
    SELECT CASE (config%preconditioner)
    CASE (NONE_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, identity_sp, &
           exchange, tol, maxiter)
    CASE (JACOBI_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, jacobi_sp, &
           exchange, tol, maxiter)
    CASE (ILU0_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, ilu0_sp, &
           exchange, tol, maxiter)
    CASE (SSOR_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, SSOR_sp, &
           exchange, tol, maxiter)
    CASE (ICC_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, icc_sp, &
           exchange, tol, maxiter)
    CASE (MICC_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, icc_sp, &
           exchange, tol, maxiter)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "] does not exist!", 'solver_all_multi.f90', 176)
    END SELECT
  CASE DEFAULT
    CALL abort_ppm("Selected solver [" // int2str(config%solver) &
         // "]  does not exist!", 'solver_all_multi.f90', 180)
  END SELECT

END FUNCTION solve_sp

!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
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
!
! Dispatcher solve function which calls the appropiate solver/preconditioner
! according to the config settings
FUNCTION solve_dp(A, b, x, ext_x, exchange, tol_opt, maxiter_opt) RESULT(kiter)
  USE solver_internal
  USE ppm_extents, ONLY: extent
  USE ppm_base, ONLY: abort_ppm

  REAL(dp), INTENT(IN) :: b(:,:)                ! right-hand-side b
  REAL(dp), INTENT(INOUT) :: x(:,:)             ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(dp), OPTIONAL, INTENT(IN) :: tol_opt     ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
  INTEGER :: kiter, maxiter
  REAL(dp) :: tol, lambda_min, lambda_max

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
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
    tol = REAL(config%tol, dp)
  ENDIF

  ! Make sure stencil is set up
  IF (.NOT. stencil_defined(0.0_precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         'solver_all_multi.f90', 91)
  ENDIF

  ! Make sure preconditioner is initialized
  IF (.NOT. precond_prepared(config%preconditioner, 0.0_precs)) THEN
    SELECT CASE (config%preconditioner)
    CASE (JACOBI_PRECOND)
      CALL prep_jacobi(0.0_precs)
    CASE (ILU0_PRECOND)
      CALL prep_ilu0(0.0_precs)
    CASE (SSOR_PRECOND)
      CALL abort_ppm("No SSOR parameter provided!", 'solver_all_multi.f90', 102)
    CASE (ICC_PRECOND)
      CALL prep_icc(config%icc_param, 0.0_precs)
    CASE (MICC_PRECOND)
      CALL prep_micc(config%icc_param, 0.0_precs)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "]  does not exist!", 'solver_all_multi.f90', 110)
    END SELECT
  ENDIF

  ! Call right preconditioner/solver combination
  SELECT CASE (config%solver)
  CASE (CG_SOLVER)
    SELECT CASE (config%preconditioner)
    CASE (NONE_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, identity_dp, &
           exchange, tol, maxiter)
    CASE (JACOBI_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, jacobi_dp, &
           exchange, tol, maxiter)
    CASE (ILU0_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, ilu0_dp, &
           exchange, tol, maxiter)
    CASE (SSOR_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, SSOR_dp, &
           exchange, tol, maxiter)
    CASE (ICC_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, icc_dp, &
           exchange, tol, maxiter)
    CASE (MICC_PRECOND)
      kiter = precond_cg_method(A, b, x, ext_x, icc_dp, &
           exchange, tol, maxiter)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "] does not exist!", 'solver_all_multi.f90', 139)
    END SELECT
  CASE (CHEBYSHEV_SOLVER)
    lambda_min = REAL(config%lambda_min, dp)
    lambda_max = REAL(config%lambda_max, dp)
    IF (ieee_is_nan(lambda_min) .OR. ieee_is_nan(lambda_max)) THEN
      CALL abort_ppm("Lambda_min and lambda_max not provided!", 'solver_all_multi.f90', &
           146)
    ENDIF
    SELECT CASE (config%preconditioner)
    CASE (NONE_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, identity_dp, &
           exchange, tol, maxiter)
    CASE (JACOBI_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, jacobi_dp, &
           exchange, tol, maxiter)
    CASE (ILU0_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, ilu0_dp, &
           exchange, tol, maxiter)
    CASE (SSOR_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, SSOR_dp, &
           exchange, tol, maxiter)
    CASE (ICC_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, icc_dp, &
           exchange, tol, maxiter)
    CASE (MICC_PRECOND)
      kiter = precond_chebyshev_method(A, b, x, ext_x, &
           lambda_min, lambda_max, icc_dp, &
           exchange, tol, maxiter)
    CASE DEFAULT
      CALL abort_ppm("Selected preconditioner [" &
           // int2str(config%preconditioner) &
           // "] does not exist!", 'solver_all_multi.f90', 176)
    END SELECT
  CASE DEFAULT
    CALL abort_ppm("Selected solver [" // int2str(config%solver) &
         // "]  does not exist!", 'solver_all_multi.f90', 180)
  END SELECT

END FUNCTION solve_dp

!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:

END MODULE solver_all
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
