!>
!! @file solver_config.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!
! Keywords: scales ppm solver configuration
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
!> configuration for stencil and parameters
#include "fc_feature_defs.inc"
MODULE solver_config
  USE solver_internal
  USE ppm_std_type_kinds, ONLY: wp, sp, dp
  USE ppm_extents, ONLY: extent, extent_start, extent_end
  USE ppm_base, ONLY: abort_ppm
  USE ppm_math_extensions, ONLY: assign_nan

  IMPLICIT NONE

  PRIVATE

  !> @brief determines if stencil is defined
  !! @param[in] prec float prototype to determine the precision
  !! @returns boolean
  INTERFACE stencil_defined
    MODULE PROCEDURE stencil_defined_sp
    MODULE PROCEDURE stencil_defined_dp
  END INTERFACE

  !> @brief set the stencil, i.e. matrix A
  !! @param[in] zonal 2-d array target
  !! @param[in] meridional 2-d array target
  !! @param[in] central 2-d array target
  !! @param[in] ext extend list of length 3 of the arrays above
  INTERFACE set_stencil
    MODULE PROCEDURE set_stencil_sp
    MODULE PROCEDURE set_stencil_dp
  END INTERFACE

  !> @brief apply stencil, i.e. matrix vector multiplication
  !! @param[in] x input 2d array
  !! @param[out] y ouput 2d array
  INTERFACE apply_stencil
    MODULE PROCEDURE apply_stencil_sp
    MODULE PROCEDURE apply_stencil_dp
  END INTERFACE

  !> @brief apply shifted stencil, i.e. shifted matrix vector multiplication
  !! @param[in] x input 2d array
  !! @param[out] y ouput 2d array
  INTERFACE apply_shifted_stencil
    MODULE PROCEDURE apply_shifted_stencil_sp
    MODULE PROCEDURE apply_shifted_stencil_dp
  END INTERFACE

  !-------------------------------------------------------------------
  ! Module Function & Procedure Interface
  !-------------------------------------------------------------------
  PUBLIC :: init_solver, set_stencil, get_ssor_param, set_ssor_param &
       , set_lambda_min_max, stencil_defined, apply_stencil &
       , apply_shifted_stencil, get_lambda_min_max, apply_stencil_sp &
       , apply_stencil_dp, apply_shifted_stencil_sp &
       , apply_shifted_stencil_dp

  CHARACTER(len=*), PARAMETER :: filename = 'solver_config.f90'
CONTAINS

  !> @brief initialize solver with a given configuration or default
  !! configuration
  !! @param[in] solver_config the @link solver_internal::solver_config_type
  !! solver config type @endlink holding all parameters
  SUBROUTINE init_solver(solver_config)
    USE solver_public, ONLY: CG_SOLVER, ICC_PRECOND

    TYPE(solver_config_type), INTENT(IN), OPTIONAL :: solver_config

    IF ( PRESENT(solver_config) ) THEN
      config = solver_config
    ELSE
      ! Set defaults
      config%solver = CG_SOLVER
      config%preconditioner = ICC_PRECOND
      config%do_exchange = .TRUE.
      config%tol = 1.0e-13_wp
      config%maxiter = 1000
      config%schwarz_tol = 1.0e-13_wp
      config%schwarz_maxiter = 1000
      config%schwarz_relax = 1.0_wp
      config%iterref_tol = 1.0e-13_wp
      config%iterref_maxiter = 1000
      CALL assign_nan(config%lambda_min) ! NaN will be treated as uninitialized
      CALL assign_nan(config%lambda_max) ! NaN will be treated as uninitialized
      CALL assign_nan(config%ssor_param) ! NaN will be treated as uninitialized
      config%icc_param = 4
      config%lambda_maxiter = 10000
      config%lambda_tol = 1.0e-13_wp
      config%cheby_miniter = 50
      config%cheby_checkrate = 5
      config%schwarz_miniter = 20
      config%schwarz_checkrate = 5
    ENDIF

  END SUBROUTINE init_solver

  ! Set SSOR parameter in configuration
  SUBROUTINE set_ssor_param(ssor_param)
    REAL(wp), INTENT(IN) :: ssor_param

    IF ( 0.0_wp < ssor_param .AND. ssor_param < 2.0_wp ) THEN
      config%ssor_param = ssor_param
    ELSE
      CALL abort_ppm("SSOR parameter should be in (0,2)", filename, __LINE__)
    ENDIF

  END SUBROUTINE set_ssor_param

  ! Get SSOR parameter from configuration
  FUNCTION get_ssor_param() RESULT (ssor_param)
    REAL(wp) :: ssor_param

    ssor_param = config%ssor_param

  END FUNCTION get_ssor_param

  ! Set eigenvalues lambda_min, lambda_max in configuration
  SUBROUTINE set_lambda_min_max(lambda_min, lambda_max)
    REAL(wp), INTENT(IN) :: lambda_min, lambda_max

    IF ( 0.0_wp < lambda_min .AND. lambda_min < lambda_max ) THEN
      config%lambda_min = lambda_min
      config%lambda_max = lambda_max
    ELSE
      CALL abort_ppm("Condition 0 < lambda_min < lambda_max violated!", &
           filename, __LINE__)
    ENDIF

  END SUBROUTINE set_lambda_min_max

  ! Get eigenvalues lambda_min, lambda_max from configuration
  FUNCTION get_lambda_min_max() RESULT(lambdas)
    REAL(wp) :: lambdas(2)

    lambdas(1) = config%lambda_min
    lambdas(2) = config%lambda_max

  END FUNCTION get_lambda_min_max

#define PREC sp
#define STENCIL stencil_sp
#define STENCIL_DEFINED stencil_defined_sp
#define SET_STENCIL set_stencil_sp
#define APPLY_STENCIL apply_stencil_sp
#define APPLY_SHIFTED_STENCIL apply_shifted_stencil_sp
#include "solver_config_multi.f90"
#undef PREC
#undef STENCIL
#undef STENCIL_DEFINED
#undef SET_STENCIL
#undef APPLY_STENCIL
#undef APPLY_SHIFTED_STENCIL
#define PREC dp
#define STENCIL stencil_dp
#define STENCIL_DEFINED stencil_defined_dp
#define SET_STENCIL set_stencil_dp
#define APPLY_STENCIL apply_stencil_dp
#define APPLY_SHIFTED_STENCIL apply_shifted_stencil_dp
#include "solver_config_multi.f90"
#undef PREC
#undef STENCIL
#undef STENCIL_DEFINED
#undef SET_STENCIL
#undef APPLY_STENCIL
#undef APPLY_SHIFTED_STENCIL

END MODULE solver_config
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
