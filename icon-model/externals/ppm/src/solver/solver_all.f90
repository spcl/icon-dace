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
#include "fc_feature_defs.inc"
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
#ifdef BLASVES
    USE ppm_extents, ONLY: extent_size
#endif
    REAL(dp), INTENT(IN) :: b(:,:)                  ! right-hand-side b
    REAL(dp), INTENT(INOUT) :: x(:,:)               ! startvalue and result
    TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
    REAL(dp), OPTIONAL, INTENT(IN) :: tol_opt       ! tolerance for residual
    INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
    INTEGER :: jb, ib, il, jl, ie, je, totaliter, kiter, kkiter, maxiter
#ifdef BLASVES
    INTEGER :: inner_size
#endif
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
#ifdef BLASVES
    inner_size = extent_size(ext_x)
#endif

    ALLOCATE(r_dp(ie,je), r_sp(ie,je), c_dp(ie,je), c_sp(ie,je), Ax(ie,je))

    ! Initialise residual and direction
    CALL A_dp(x, Ax)
#ifndef BLASVES
    r_dp(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)
#else
    CALL DVES(inner_size,b(ib:il,jb:jl),1,Ax(ib:il,jb:jl),1,r_dp(ib:il,jb:jl),1)
#endif
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
#ifndef BLASVES
      r_dp(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)
#else
      CALL DVES(inner_size,b(ib:il,jb:jl),1,Ax(ib:il,jb:jl),1,&
           r_dp(ib:il,jb:jl),1)
#endif
      delta = arr_dotproduct(r_dp(ib:il,jb:jl))
      relres = delta*delta0
      kiter = kiter + 1
      totaliter = totaliter + kkiter

    ENDDO

    DEALLOCATE(r_dp, r_sp, c_dp, c_sp, Ax)

  END FUNCTION iter_refinement

#define PREC sp
#define SOLVE solve_sp
#define IDENTITY identity_sp
#define JACOBI jacobi_sp
#define ILU0 ilu0_sp
#define SSOR SSOR_sp
#define ICC icc_sp
#include "solver_all_multi.f90"
#undef PREC
#undef SOLVE
#undef IDENTITY
#undef JACOBI
#undef ILU0
#undef SSOR
#undef ICC
#define PREC dp
#define SOLVE solve_dp
#define IDENTITY identity_dp
#define JACOBI jacobi_dp
#define ILU0 ilu0_dp
#define SSOR SSOR_dp
#define ICC icc_dp
#include "solver_all_multi.f90"
#undef PREC
#undef SOLVE
#undef IDENTITY
#undef JACOBI
#undef ILU0
#undef SSOR
#undef ICC

END MODULE solver_all
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
