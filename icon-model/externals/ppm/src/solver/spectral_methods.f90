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
#include "fc_feature_defs.inc"
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

#define PREC sp
#define STENCIL stencil_sp
#define POWER_METHOD power_method_sp
#define CALC_LAMBDA_MAX calc_lambda_max_sp
#define CALC_LAMBDA_MIN calc_lambda_min_sp
#define CALC_LAMBDA_MIN_MAX calc_lambda_min_max_sp
#define CALC_OPTIMAL_SOR_PARAM calc_optimal_sor_param_sp
#define JACOBI_ITER jacobi_iter_sp
#include "spectral_methods_multi.f90"
#undef PREC
#undef STENCIL
#undef POWER_METHOD
#undef CALC_LAMBDA_MAX
#undef CALC_LAMBDA_MIN
#undef CALC_LAMBDA_MIN_MAX
#undef CALC_OPTIMAL_SOR_PARAM
#undef JACOBI_ITER
#define PREC dp
#define STENCIL stencil_dp
#define POWER_METHOD power_method_dp
#define CALC_LAMBDA_MAX calc_lambda_max_dp
#define CALC_LAMBDA_MIN calc_lambda_min_dp
#define CALC_LAMBDA_MIN_MAX calc_lambda_min_max_dp
#define CALC_OPTIMAL_SOR_PARAM calc_optimal_sor_param_dp
#define JACOBI_ITER jacobi_iter_dp
#include "spectral_methods_multi.f90"
#undef PREC
#undef STENCIL
#undef POWER_METHOD
#undef CALC_LAMBDA_MAX
#undef CALC_LAMBDA_MIN
#undef CALC_LAMBDA_MIN_MAX
#undef CALC_OPTIMAL_SOR_PARAM
#undef JACOBI_ITER

END MODULE spectral_methods
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
