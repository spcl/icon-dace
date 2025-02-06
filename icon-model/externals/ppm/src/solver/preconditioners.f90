!>
!! @file preconditioners.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!
! Version: 1.0
! Keywords: scales ppm solver preconditioners ilu0 icc jacobi ssor
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
!> @brief preconditioners for symmetric 5-point stencil system
!! @details all functions/subroutines have single/double precision variants
#include "fc_feature_defs.inc"
MODULE preconditioners
  USE solver_public, ONLY: micc_precond, none_precond, jacobi_precond, &
       ilu0_precond, ssor_precond, icc_precond
  USE solver_internal, ONLY: stencil_sp, stencil_dp, config
  USE ieee_arithmetic, ONLY: ieee_is_nan
  USE ppm_std_type_kinds, ONLY: wp, sp, dp
  USE ppm_base, ONLY: abort_ppm
  USE ppm_extents, ONLY: extent_start, extent_end, extent_size
  USE solver_config, ONLY: get_ssor_param, apply_stencil

  IMPLICIT NONE

  PRIVATE

  !-------------------------------------------------------------------
  ! Internal State Variables
  !-------------------------------------------------------------------
  ! Auxiliary variable for jacobi preconditioner
  REAL(sp), ALLOCATABLE :: jacobi_diag_sp(:,:)
  REAL(dp), ALLOCATABLE :: jacobi_diag_dp(:,:)
  ! Auxiliary variables for ilu0 preconditioner
  REAL(sp), ALLOCATABLE :: ilu0_diag_sp(:,:)
  REAL(dp), ALLOCATABLE :: ilu0_diag_dp(:,:)
  ! Auxiliary variables for icc(p) preconditioner
  REAL(sp), ALLOCATABLE :: ICC_C_sp(:,:), ICC_W_sp(:,:,:), ICC_S_sp(:,:,:)
  REAL(dp), ALLOCATABLE :: ICC_C_dp(:,:), ICC_W_dp(:,:,:), ICC_S_dp(:,:,:)

  !-------------------------------------------------------------------
  ! Module Function & Procedure Interface
  !-------------------------------------------------------------------

  !> @brief identity function
  !! @details use this to obtain the non-preconditioned CG-method
  !! @param[in,out] r current residual, 2d array
  INTERFACE identity
    MODULE PROCEDURE identity_sp
    MODULE PROCEDURE identity_dp
  END INTERFACE

  !> @brief prepare Jacobi preconditioner
  !! @param[in] prec float prototype to determine the precision
  INTERFACE prep_jacobi
    MODULE PROCEDURE prep_jacobi_sp
    MODULE PROCEDURE prep_jacobi_dp
  END INTERFACE

  !> @brief Jacobi preconditioner
  !! @param[in,out] r current residual, 2d array
  INTERFACE jacobi
    MODULE PROCEDURE jacobi_sp
    MODULE PROCEDURE jacobi_dp
  END INTERFACE

  !> @brief prepare ILU(0) preconditioner
  !! @param[in] prec float prototype to determine the precision
  INTERFACE prep_ilu0
    MODULE PROCEDURE prep_ilu0_sp
    MODULE PROCEDURE prep_ilu0_dp
  END INTERFACE

  !> @brief incomplete LU-decomposition with fill-in 0
  !! @param[in,out] r current residual, 2d array
  INTERFACE ilu0
    MODULE PROCEDURE ilu0_sp
    MODULE PROCEDURE ilu0_dp
  END INTERFACE

  !> @brief symmetric, successive over-relaxation preconditioner
  !! @param[in,out] r current residual, 2d array
  INTERFACE ssor
    MODULE PROCEDURE ssor_sp
    MODULE PROCEDURE ssor_dp
  END INTERFACE

  !> @brief prepare ICC(p) preconditioner
  !! @param[in] p level of fill-in
  !! @param[in] prec float prototype to determine the precision
  !! @param[in] modify boolean, use modified ICC [default: FALSE]
  INTERFACE prep_icc
    MODULE PROCEDURE prep_icc_sp
    MODULE PROCEDURE prep_icc_dp
  END INTERFACE

  !> @brief (modified) incomplete Cholesky preconditioner with fill-in p
  !! @param[in,out] r current residual, 2d array
  INTERFACE icc
    MODULE PROCEDURE icc_sp
    MODULE PROCEDURE icc_dp
  END INTERFACE

  !> @brief prepare modified ICC(p) preconditioner
  !! @param[in] p level of fill-in
  !! @param[in] prec float prototype to determine the precision
  INTERFACE prep_micc
    MODULE PROCEDURE prep_micc_sp
    MODULE PROCEDURE prep_micc_dp
  END INTERFACE

  !> @brief determines if given preconditioner is prepared
  !! @param[in] precond parameter as defined in @link solver_public @endlink
  !! @param[in] prec float prototype to determine the precision
  !! @return boolean
  INTERFACE precond_prepared
    MODULE PROCEDURE precond_prepared_sp
    MODULE PROCEDURE precond_prepared_dp
  END INTERFACE

  !> @brief the stencil of the Jacobi preconditioned original stencil
  !! @details let M be the preconditioner and A the original matrix, then
  !! this functions is M^-1*A. This is used to calculate the eigenvalues of
  !! M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE jacobi_precond_stencil
    MODULE PROCEDURE jacobi_precond_stencil_sp
    MODULE PROCEDURE jacobi_precond_stencil_dp
  END INTERFACE

  !> @brief the shifted stencil of the Jacobi preconditioned original stencil
  !! @details let M be the preconditioner, I the unity matrix and A the original
  !! matrix, then this functions is (s*I - M^-1*A) with the scalar shift s.
  !! This is used to calculate the smallest absolute eigenvalues of M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE jacobi_precond_shifted_stencil
    MODULE PROCEDURE jacobi_precond_shifted_stencil_sp
    MODULE PROCEDURE jacobi_precond_shifted_stencil_dp
  END INTERFACE

  !> @brief the stencil of the ILU(0) preconditioned original stencil
  !! @details let M be the preconditioner and A the original matrix, then
  !! this functions is M^-1*A. This is used to calculate the eigenvalues of
  !! M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE ilu0_precond_stencil
    MODULE PROCEDURE ilu0_precond_stencil_sp
    MODULE PROCEDURE ilu0_precond_stencil_dp
  END INTERFACE

  !> @brief the shifted stencil of the ILU(0) preconditioned original stencil
  !! @details let M be the preconditioner, I the unity matrix and A the original
  !! matrix, then this functions is (s*I - M^-1*A) with the scalar shift s.
  !! This is used to calculate the smallest absolute eigenvalues of M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE ilu0_precond_shifted_stencil
    MODULE PROCEDURE ilu0_precond_shifted_stencil_sp
    MODULE PROCEDURE ilu0_precond_shifted_stencil_dp
  END INTERFACE

  !> @brief the stencil of the SSOR preconditioned original stencil
  !! @details let M be the preconditioner and A the original matrix, then
  !! this functions is M^-1*A. This is used to calculate the eigenvalues of
  !! M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE ssor_precond_stencil
    MODULE PROCEDURE ssor_precond_stencil_sp
    MODULE PROCEDURE ssor_precond_stencil_dp
  END INTERFACE

  !> @brief the shifted stencil of the SSOR preconditioned original stencil
  !! @details let M be the preconditioner, I the unity matrix and A the original
  !! matrix, then this functions is (s*I - M^-1*A) with the scalar shift s.
  !! This is used to calculate the smallest absolute eigenvalues of M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE ssor_precond_shifted_stencil
    MODULE PROCEDURE ssor_precond_shifted_stencil_sp
    MODULE PROCEDURE ssor_precond_shifted_stencil_dp
  END INTERFACE

  !> @brief the stencil of the (modified) ICC preconditioned original stencil
  !! @details let M be the preconditioner and A the original matrix, then
  !! this functions is M^-1*A. This is used to calculate the eigenvalues of
  !! M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE icc_precond_stencil
    MODULE PROCEDURE icc_precond_stencil_sp
    MODULE PROCEDURE icc_precond_stencil_dp
  END INTERFACE

  !> @brief the shifted stencil of the (modified) ICC preconditioned original
  !! stencil
  !! @details let M be the preconditioner, I the unity matrix and A the original
  !! matrix, then this functions is (s*I - M^-1*A) with the scalar shift s.
  !! This is used to calculate the smallest absolute eigenvalues of M^-1*A.
  !! @param[in] x 2d-field
  !! @param[out] y 2d-field
  INTERFACE icc_precond_shifted_stencil
    MODULE PROCEDURE icc_precond_shifted_stencil_sp
    MODULE PROCEDURE icc_precond_shifted_stencil_dp
  END INTERFACE

  PUBLIC :: NONE_PRECOND, JACOBI_PRECOND, ILU0_PRECOND, SSOR_PRECOND &
       , ICC_PRECOND, MICC_PRECOND
  PUBLIC :: identity, identity_sp, identity_dp, prep_jacobi, jacobi &
       , jacobi_sp, jacobi_dp, prep_ilu0, ilu0, ilu0_sp, ilu0_dp, ssor &
       , ssor_sp, ssor_dp, prep_icc, icc, icc_sp, icc_dp &
       , precond_prepared, jacobi_precond_stencil, ilu0_precond_stencil &
       , ssor_precond_stencil, icc_precond_stencil &
       , jacobi_precond_shifted_stencil, ilu0_precond_shifted_stencil &
       , ssor_precond_shifted_stencil, icc_precond_shifted_stencil &
       , prep_micc &
       , jacobi_precond_stencil_sp, jacobi_precond_stencil_dp &
       , jacobi_precond_shifted_stencil_sp &
       , jacobi_precond_shifted_stencil_dp &
       , ilu0_precond_stencil_sp, ilu0_precond_stencil_dp &
       , ilu0_precond_shifted_stencil_sp &
       , ilu0_precond_shifted_stencil_dp &
       , ssor_precond_stencil_sp, ssor_precond_stencil_dp &
       , ssor_precond_shifted_stencil_sp &
       , ssor_precond_shifted_stencil_dp &
       , icc_precond_stencil_sp, icc_precond_stencil_dp &
       , icc_precond_shifted_stencil_sp &
       , icc_precond_shifted_stencil_dp

CONTAINS

#define PREC sp
#define STENCIL stencil_sp
#define IDENTITY identity_sp
#define JACOBI_DIAG jacobi_diag_sp
#define PREP_JACOBI prep_jacobi_sp
#define JACOBI jacobi_sp
#define ILU0_DIAG ilu0_diag_sp
#define PREP_ILU0 prep_ilu0_sp
#define ILU0 ilu0_sp
#define SSOR SSOR_sp
#define ICC_C ICC_C_sp
#define ICC_W ICC_W_sp
#define ICC_S ICC_S_sp
#define PREP_ICC prep_icc_sp
#define ICC icc_sp
#define PREP_MICC prep_micc_sp
#define PRECOND_PREPARED precond_prepared_sp
#define JACOBI_PRECOND_STENCIL jacobi_precond_stencil_sp
#define JACOBI_PRECOND_SHIFTED_STENCIL jacobi_precond_shifted_stencil_sp
#define ILU0_PRECOND_STENCIL ilu0_precond_stencil_sp
#define ILU0_PRECOND_SHIFTED_STENCIL ilu0_precond_shifted_stencil_sp
#define SSOR_PRECOND_STENCIL ssor_precond_stencil_sp
#define SSOR_PRECOND_SHIFTED_STENCIL ssor_precond_shifted_stencil_sp
#define ICC_PRECOND_STENCIL icc_precond_stencil_sp
#define ICC_PRECOND_SHIFTED_STENCIL icc_precond_shifted_stencil_sp
#include "preconditioners_multi.f90"
#undef PREC
#undef STENCIL
#undef IDENTITY
#undef JACOBI_DIAG
#undef PREP_JACOBI
#undef JACOBI
#undef ILU0_DIAG
#undef PREP_ILU0
#undef ILU0
#undef SSOR
#undef ICC_C
#undef ICC_W
#undef ICC_S
#undef PREP_ICC
#undef ICC
#undef PREP_MICC
#undef PRECOND_PREPARED
#undef JACOBI_PRECOND_STENCIL
#undef JACOBI_PRECOND_SHIFTED_STENCIL
#undef ILU0_PRECOND_STENCIL
#undef ILU0_PRECOND_SHIFTED_STENCIL
#undef SSOR_PRECOND_STENCIL
#undef SSOR_PRECOND_SHIFTED_STENCIL
#undef ICC_PRECOND_STENCIL
#undef ICC_PRECOND_SHIFTED_STENCIL
#define PREC dp
#define STENCIL stencil_dp
#define IDENTITY identity_dp
#define JACOBI_DIAG jacobi_diag_dp
#define PREP_JACOBI prep_jacobi_dp
#define JACOBI jacobi_dp
#define ILU0_DIAG ilu0_diag_dp
#define PREP_ILU0 prep_ilu0_dp
#define ILU0 ilu0_dp
#define SSOR SSOR_dp
#define ICC_C ICC_C_dp
#define ICC_W ICC_W_dp
#define ICC_S ICC_S_dp
#define PREP_ICC prep_icc_dp
#define ICC icc_dp
#define PREP_MICC prep_micc_dp
#define PRECOND_PREPARED precond_prepared_dp
#define JACOBI_PRECOND_STENCIL jacobi_precond_stencil_dp
#define JACOBI_PRECOND_SHIFTED_STENCIL jacobi_precond_shifted_stencil_dp
#define ILU0_PRECOND_STENCIL ilu0_precond_stencil_dp
#define ILU0_PRECOND_SHIFTED_STENCIL ilu0_precond_shifted_stencil_dp
#define SSOR_PRECOND_STENCIL ssor_precond_stencil_dp
#define SSOR_PRECOND_SHIFTED_STENCIL ssor_precond_shifted_stencil_dp
#define ICC_PRECOND_STENCIL icc_precond_stencil_dp
#define ICC_PRECOND_SHIFTED_STENCIL icc_precond_shifted_stencil_dp
#include "preconditioners_multi.f90"
#undef PREC
#undef STENCIL
#undef IDENTITY
#undef JACOBI_DIAG
#undef PREP_JACOBI
#undef JACOBI
#undef ILU0_DIAG
#undef PREP_ILU0
#undef ILU0
#undef SSOR
#undef ICC_C
#undef ICC_W
#undef ICC_S
#undef PREP_ICC
#undef ICC
#undef PREP_MICC
#undef PRECOND_PREPARED
#undef JACOBI_PRECOND_STENCIL
#undef JACOBI_PRECOND_SHIFTED_STENCIL
#undef ILU0_PRECOND_STENCIL
#undef ILU0_PRECOND_SHIFTED_STENCIL
#undef SSOR_PRECOND_STENCIL
#undef SSOR_PRECOND_SHIFTED_STENCIL
#undef ICC_PRECOND_STENCIL
#undef ICC_PRECOND_SHIFTED_STENCIL

END MODULE preconditioners
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
