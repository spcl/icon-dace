!>
!! @file solvers.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!! @author Florian Wilhelm <Florian.Wilhelm@kit.edu>
!
! Keywords: scales ppm solver cg chebyshev schwarz
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
!> provides solvers i.e. preconditioned CG, Chebyshev, Schwarz
#include "fc_feature_defs.inc"
MODULE solvers
  USE solver_public
  USE solver_internal
  USE linear_algebra
  USE ppm_std_type_kinds, ONLY: sp, dp
  USE ppm_extents, ONLY: extent_start, extent_end, extent
#if defined BLASVES || defined BLASAXPY || defined BLASSCAL
  USE ppm_extents, ONLY: extent_size
#endif
  USE preconditioners, ONLY: identity_sp, identity_dp

  IMPLICIT NONE

  PRIVATE

  !-------------------------------------------------------------------
  ! Module Function & Procedure Interface
  !-------------------------------------------------------------------

  !> @brief preconditioned CG method
  !! @details variants in single and double precision
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] b right hand side
  !! @param[in,out] x start vector, holds solution in the end
  !! @param[in] ext_x extent of x
  !! @param[in] precond preconditioner function from @link preconditioners
  !! @endlink
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return number of iterations
  INTERFACE precond_cg_method
    MODULE PROCEDURE precond_cg_method_sp
    MODULE PROCEDURE precond_cg_method_dp
  END INTERFACE

  !> @brief CG method
  !! @details variants in single and double precision
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
  INTERFACE cg_method
    MODULE PROCEDURE cg_method_sp
    MODULE PROCEDURE cg_method_dp
  END INTERFACE

  !> @brief preconditioned Chebyshev iteration
  !! @details variants in single and double precision
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] b right hand side
  !! @param[in,out] x start vector, holds solution in the end
  !! @param[in] ext_x extent of x
  !! @param[in] lambda_min smallest absolute eigenvalue
  !! @param[in] lambda_max largest absolute eigenvalue
  !! @param[in] precond preconditioner function from @link preconditioners
  !! @endlink
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return number of iterations
  INTERFACE precond_chebyshev_method
    MODULE PROCEDURE precond_chebyshev_method_sp
    MODULE PROCEDURE precond_chebyshev_method_dp
  END INTERFACE

  !> @brief Chebyshev iteration
  !! @details variants in single and double precision
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] b right hand side
  !! @param[in,out] x start vector, holds solution in the end
  !! @param[in] ext_x extent of x
  !! @param[in] lambda_min smallest absolute eigenvalue
  !! @param[in] lambda_max largest absolute eigenvalue
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @return number of iterations
  INTERFACE chebyshev_method
    MODULE PROCEDURE chebyshev_method_sp
    MODULE PROCEDURE chebyshev_method_dp
  END INTERFACE

  !> @brief Additive Schwarz method
  !! @details variants in single and double precision
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] b right hand side
  !! @param[in,out] x start vector, holds solution in the end
  !! @param[in] ext_x extent of x
  !! @param[in] lambda_min smallest absolute eigenvalue
  !! @param[in] lambda_max largest absolute eigenvalue
  !! @param[in] local_solver solver function like CG method
  !! @param[in] exchange @link solver_internal::exchangeop
  !! boundary exchange @endlink function
  !! @param[in] tol_opt optional stopping tolerance for relative residual
  !! @param[in] maxiter_opt optional maximum number of iterations
  !! @param[in] omega_opt optional relaxation parameter
  !! @return number of iterations
  INTERFACE schwarz_method
    MODULE PROCEDURE schwarz_method_sp
    MODULE PROCEDURE schwarz_method_dp
  END INTERFACE

  PUBLIC :: CG_SOLVER, CHEBYSHEV_SOLVER
  PUBLIC :: cg_method, precond_cg_method, chebyshev_method, schwarz_method, &
       precond_chebyshev_method

CONTAINS

#define PREC sp
#define PAXPY SAXPY
#define PVES SVES
#define PYAX SYAX
#define PSCAL SSCAL
#define IDENTITY identity_sp
#define PRECOND_CG_METHOD precond_cg_method_sp
#define CG_METHOD cg_method_sp
#define SCHWARZ_METHOD schwarz_method_sp
#define CHEBYSHEV_METHOD chebyshev_method_sp
#define PRECOND_CHEBYSHEV_METHOD precond_chebyshev_method_sp
#include "solvers_multi.f90"
#undef PREC
#undef PAXPY
#undef PVES
#undef PYAX
#undef PSCAL
#undef IDENTITY
#undef PRECOND_CG_METHOD
#undef CG_METHOD
#undef SCHWARZ_METHOD
#undef CHEBYSHEV_METHOD
#undef PRECOND_CHEBYSHEV_METHOD
#define PREC dp
#define PAXPY DAXPY
#define PVES DVES
#define PYAX DYAX
#define PSCAL DSCAL
#define IDENTITY identity_dp
#define PRECOND_CG_METHOD precond_cg_method_dp
#define CG_METHOD cg_method_dp
#define SCHWARZ_METHOD schwarz_method_dp
#define CHEBYSHEV_METHOD chebyshev_method_dp
#define PRECOND_CHEBYSHEV_METHOD precond_chebyshev_method_dp
#include "solvers_multi.f90"
#undef PREC
#undef PAXPY
#undef PVES
#undef PYAX
#undef PSCAL
#undef IDENTITY
#undef PRECOND_CG_METHOD
#undef CG_METHOD
#undef SCHWARZ_METHOD
#undef CHEBYSHEV_METHOD
#undef PRECOND_CHEBYSHEV_METHOD
END MODULE solvers
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
