!>
!! @file linear_algebra.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!! @version 1.0
!
! Keywords: scales ppm solver linear algebra residual
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
!> basic tools from linear algebra
#include "fc_feature_defs.inc"
MODULE linear_algebra
  USE solver_internal, ONLY: config
  USE ppm_std_type_kinds, ONLY: sp, dp
  USE ppm_extents, ONLY: extent, extent_start, extent_end
#if defined USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif

  PRIVATE

  !-------------------------------------------------------------------
  ! Module Function & Procedure Interface
  !-------------------------------------------------------------------

  !> @brief Dot product for dim 1 and 2 arrays in single/double precision.
  !! @details If only one array is povided perform dot product x^T * x.
  !! @param[in] x array
  !! @param[in] y optional second array
  !! @param[in] global_opt boolean, perform operation globally if true
  !! @return dot product
  INTERFACE arr_dotproduct
    MODULE PROCEDURE arr_dotproduct1_sp
    MODULE PROCEDURE arr_dotproduct1_dp
    MODULE PROCEDURE arr_dotproduct2_sp
    MODULE PROCEDURE arr_dotproduct2_dp
  END INTERFACE

  !> @brief Calculate absolute residual in single/double precision
  !! @details Blballala
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] b right hand side
  !! @param[in] x approximate solution x
  !! @param[in] ext_x extent of x
  !! @param[in] global_opt boolean, perform operation globally if true
  !! @return global/local absolute residual
  INTERFACE calc_abs_res
    MODULE PROCEDURE calc_abs_res_sp
    MODULE PROCEDURE calc_abs_res_dp
  END INTERFACE

  !> @brief Calculate relative residual in single/double precision for Ax=b
  !! @param[in] A @link solver_internal::linop linear operator@endlink like
  !! matrix multiplication
  !! @param[in] b right hand side
  !! @param[in] x approximate solution x
  !! @param[in] ext_x extent of x
  !! @param[in] global_opt boolean, perform operation globally if true
  !! @return global/local absolute residual
  INTERFACE calc_rel_res
    MODULE PROCEDURE calc_rel_res_sp
    MODULE PROCEDURE calc_rel_res_dp
  END INTERFACE

  !> @brief Dot product for dim 1 and 2 arrays in single/double precision.
  !! @details If only one array is povided perform dot product x^T * x.
  !! @param[in] x array
  !! @param[in] y optional second array
  !! @param[in] global_opt boolean, perform operation globally if true
  !! @return dot product
  INTERFACE global_sum
    MODULE PROCEDURE global_sum_sp
    MODULE PROCEDURE global_sum_dp
  END INTERFACE

  !> @brief The 2-norm of a 2d array
  !! @param[in] x 2d array
  !! @param[in] global_opt boolean, perform operation globally if true
  !! @return dot product
  INTERFACE arr_norm_2
    MODULE PROCEDURE arr_norm_2_sp
    MODULE PROCEDURE arr_norm_2_dp
  END INTERFACE

  PUBLIC :: calc_abs_res, calc_rel_res, arr_dotproduct, global_sum, arr_norm_2

#if defined __GNUC__ && __GNUC__ > 4 \
  && defined USE_MPI && ! defined USE_MPI_MOD
  INTERFACE
    SUBROUTINE mpi_allreduce(sendbuf, recvbuf, count, datatype, op, comm, &
         ierror)
      INTEGER, INTENT(in) :: count, datatype, op, comm
      INTEGER, INTENT(out) :: ierror
      INTEGER, INTENT(in) :: sendbuf(*)
      INTEGER, INTENT(inout) :: recvbuf(*)
!DEC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
!GCC$ ATTRIBUTES NO_ARG_CHECK :: sendbuf, recvbuf
!$PRAGMA IGNORE_TKR sendbuf, recvbuf
!DIR$ IGNORE_TKR sendbuf, recvbuf
!IBM* IGNORE_TKR sendbuf, recvbuf
    END SUBROUTINE mpi_allreduce
  END INTERFACE
#endif

CONTAINS

#define PREC sp
#define PDOT SDOT
#define PREC_MPI_DT mp_sp
#define ARR_DOTPRODCUT1 arr_dotproduct1_sp
#define ARR_DOTPRODCUT2 arr_dotproduct2_sp
#define CALC_ABS_RES calc_abs_res_sp
#define CALC_REL_RES calc_rel_res_sp
#define GLOBAL_SUM global_sum_sp
#define ARR_NORM_2 arr_norm_2_sp
#include "linear_algebra_multi.f90"
#undef PREC
#undef PDOT
#undef PREC_MPI_DT
#undef ARR_DOTPRODCUT1
#undef ARR_DOTPRODCUT2
#undef CALC_ABS_RES
#undef CALC_REL_RES
#undef GLOBAL_SUM
#undef ARR_NORM_2
#define PREC dp
#define PDOT DDOT
#define PREC_MPI_DT mp_dp
#define ARR_DOTPRODCUT1 arr_dotproduct1_dp
#define ARR_DOTPRODCUT2 arr_dotproduct2_dp
#define CALC_ABS_RES calc_abs_res_dp
#define CALC_REL_RES calc_rel_res_dp
#define GLOBAL_SUM global_sum_dp
#define ARR_NORM_2 arr_norm_2_dp
#include "linear_algebra_multi.f90"
#undef PREC
#undef PDOT
#undef PREC_MPI_DT
#undef ARR_DOTPRODCUT1
#undef ARR_DOTPRODCUT2
#undef CALC_ABS_RES
#undef CALC_REL_RES
#undef GLOBAL_SUM
#undef ARR_NORM_2

END MODULE linear_algebra
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
