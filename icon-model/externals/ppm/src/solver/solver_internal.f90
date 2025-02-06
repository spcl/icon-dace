!>
!! @file solver_internal.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!
! Keywords: scales ppm solver internal
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
!> @brief internal solver module
!! @details Avoid including/importing this module as it is supposed to be only
!! used from other solver modules.
#include "fc_feature_defs.inc"
MODULE solver_internal
  USE ppm_std_type_kinds, ONLY: wp, sp, dp
  USE ppm_extents, ONLY: extent, extent_start, extent_end
  USE ppm_base, ONLY: abort_ppm
  USE ieee_arithmetic, ONLY: is_normal => ieee_is_normal

#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif

  !> @brief the stencil type, i.e. structure in single precision
  !! @details it holds pointers to the zonal, meridional and central components
  !! as well as a list of extents for the target arrays. Additionally, a shift s
  !! is given to calculate (s*I - A) with unity matrix I and stencil matrix A.
  TYPE stencil_type_sp
    REAL(sp), DIMENSION(:,:), POINTER :: zonal => NULL(), meridional => NULL(),&
         central => NULL()
    TYPE(extent), DIMENSION(:), POINTER :: extent => NULL()
    REAL(sp) :: shift = 0.0_sp ! For shifting when determing least dominant EV
  END TYPE stencil_type_sp

  !> @brief the stencil type, i.e. structure in double precision
  !! @details it holds pointers to the zonal, meridional and central components
  !! as well as a list of extents for the target arrays. Additionally, a shift s
  !! is given to calculate (s*I - A) with unity matrix I and stencil matrix A.
  TYPE stencil_type_dp
    REAL(dp), DIMENSION(:,:), POINTER :: zonal => NULL(), meridional => NULL(),&
         central => NULL()
    TYPE(extent), DIMENSION(:), POINTER :: extent => NULL()
    REAL(dp) :: shift = 0.0_dp ! For shifting when determing least dominant EV
  END TYPE stencil_type_dp

  !> type holding all configuration parameters
  TYPE solver_config_type
    !> solver to use as defined as parameter in @link solver_public @endlink
    INTEGER :: solver
    !> which preconditioner to use as defined in @link solver_public @endlink
    INTEGER :: preconditioner
    !> perform boundary exchanges or not
    LOGICAL :: do_exchange
    !> relative residual tolerance for solver
    REAL(wp):: tol
    !> maximum number of iteration to perform
    INTEGER :: maxiter
    !> relative residual tolerance for schwarz method
    REAL(wp):: schwarz_tol
    !> maximum number of iteration to perform for Schwarz method
    INTEGER :: schwarz_maxiter
    !> relaxation parameter of Schwarz method
    REAL(wp):: schwarz_relax
    !> relative residual tolerance for iterative refinement
    REAL(wp):: iterref_tol
    !> maximum number of iterations for multi-precision iterative refinement
    INTEGER :: iterref_maxiter
    !> minimum absolute eigenvalue, take 0.0 to determine automatically
    REAL(wp):: lambda_min
    !> maximum absolute eigenvalue, take 0 to determine automatically
    REAL(wp):: lambda_max
    !> tolerance when determing lambda_min and lambda_max
    REAL(wp):: lambda_tol
    !> maximum number of iterations when determing eigenvalues
    INTEGER :: lambda_maxiter
    !> SSOR preconditioner parameter
    REAL(wp):: ssor_param
    !> level of fill-in for (modified) ICC preconditioner
    INTEGER :: icc_param
    !> minimum number of iterations for Chebyshev iteration
    INTEGER :: cheby_miniter
    !> rate of checking for break condition in Chebyshev iteration
    INTEGER :: cheby_checkrate
    !> minimum number of iterations for Schwarz method
    INTEGER :: schwarz_miniter
    !> rate of chaking for break condition in Schwarz method
    INTEGER :: schwarz_checkrate
  END TYPE solver_config_type

  !-------------------------------------------------------------------
  ! Internal State Variables
  !-------------------------------------------------------------------
  TYPE(stencil_type_sp), SAVE :: stencil_sp
  TYPE(stencil_type_dp), SAVE :: stencil_dp
  TYPE(solver_config_type), SAVE :: config

  !> @brief checks each element of a 0 to 3 dimensional array for NaN, Infty...
  !! @details variants in single and double precision
  !! @param[in] x scalar or array of maximum dimension 3
  !! @param[in] str string that is shown when non-normal number is found
  INTERFACE abort_unless_normal
    MODULE PROCEDURE abort_unless_normal0_sp
    MODULE PROCEDURE abort_unless_normal0_dp
    MODULE PROCEDURE abort_unless_normal1_sp
    MODULE PROCEDURE abort_unless_normal1_dp
    MODULE PROCEDURE abort_unless_normal2_sp
    MODULE PROCEDURE abort_unless_normal2_dp
    MODULE PROCEDURE abort_unless_normal3_sp
    MODULE PROCEDURE abort_unless_normal3_dp
  END INTERFACE

  !> @brief clears the halos of a 2d array
  !! @details variants in single and double precision
  !! @param[in] x 2d array
  !! @param[in] ext_x extent of x
  INTERFACE clear_halos
    MODULE PROCEDURE clear_halos_sp
    MODULE PROCEDURE clear_halos_dp
  END INTERFACE

  INTERFACE
    !> @brief stub function for linear operator
    !! @param[in] x input 2d array
    !! @param[out] y output 2d array
    SUBROUTINE linop(x, y)
      IMPORT :: wp
      REAL(wp), INTENT(IN)  :: x(:,:)
      REAL(wp), INTENT(OUT) :: y(:,:)
    END SUBROUTINE linop
    !> @brief stub boundary exchange operation
    !! @param[in, out] x 2d array
    !! @param[in] str optional string
    SUBROUTINE exchangeop(x, str)
      IMPORT :: wp
      REAL(wp), INTENT(INOUT) :: x(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: str
    END SUBROUTINE exchangeop
  END INTERFACE

  PUBLIC :: stencil_type_sp, stencil_type_dp, solver_config_type, stencil_sp, &
       stencil_dp, config, int2str, abort_unless_normal, clear_halos

CONTAINS

  !> @brief returns string representation of given integer i
  !! @param[in] i integer number
  !! @returns string
  PURE FUNCTION int2str(i) RESULT(cout)
    INTEGER, INTENT(IN) :: i
    CHARACTER(len=32) :: cin
    CHARACTER(len=FLOOR(LOG10(REAL(MAX(1,ABS(i)))))+1-MIN(SIGN(1,i),0)) :: cout

    WRITE(cin,"(I32)") i
    cout = ADJUSTL(cin)

  END FUNCTION int2str

  !> @brief lets rank 0 process output given strings
  !! @param[in] t1 first string
  !! @param[in] t2 optional second string
  !! @param[in] t3 optional third string
  !! @param[in] t4 optional fourth string
  !! @param[in] t5 optional fifth string
  SUBROUTINE debug(t1, t2, t3, t4, t5)
    CHARACTER(len=*), INTENT(IN) :: t1
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: t2, t3, t4, t5

    INTEGER :: my_rank
#ifdef USE_MPI
    INTEGER :: ierror
#endif

#ifdef USE_MPI
    CALL MPI_Comm_rank(mpi_comm_world, my_rank, ierror)
#else
    my_rank = 0
#endif

    IF (my_rank == 0) THEN
      IF (PRESENT(t5)) THEN
        WRITE(0,*) t1, t2, t3, t4, t5
      ELSEIF (PRESENT(t4)) THEN
        WRITE(0,*) t1, t2, t3, t4
      ELSEIF (PRESENT(t3)) THEN
        WRITE(0,*) t1, t2, t3
      ELSEIF (PRESENT(t2)) THEN
        WRITE(0,*) t1, t2
      ELSE
        WRITE(0,*) t1
      ENDIF
    ENDIF
  END SUBROUTINE debug

  ! Include SP and DP variants here
#define PREC sp
#define ABORT_UNLESS_NORMAL0 abort_unless_normal0_sp
#define ABORT_UNLESS_NORMAL1 abort_unless_normal1_sp
#define ABORT_UNLESS_NORMAL2 abort_unless_normal2_sp
#define ABORT_UNLESS_NORMAL3 abort_unless_normal3_sp
#define CLEAR_HALOS clear_halos_sp
#include "solver_internal_multi.f90"
#undef PREC
#undef ABORT_UNLESS_NORMAL0
#undef ABORT_UNLESS_NORMAL1
#undef ABORT_UNLESS_NORMAL2
#undef ABORT_UNLESS_NORMAL3
#undef CLEAR_HALOS
#define PREC dp
#define ABORT_UNLESS_NORMAL0 abort_unless_normal0_dp
#define ABORT_UNLESS_NORMAL1 abort_unless_normal1_dp
#define ABORT_UNLESS_NORMAL2 abort_unless_normal2_dp
#define ABORT_UNLESS_NORMAL3 abort_unless_normal3_dp
#define CLEAR_HALOS clear_halos_dp
#include "solver_internal_multi.f90"
#undef PREC
#undef ABORT_UNLESS_NORMAL0
#undef ABORT_UNLESS_NORMAL1
#undef ABORT_UNLESS_NORMAL2
#undef ABORT_UNLESS_NORMAL3
#undef CLEAR_HALOS

END MODULE solver_internal
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
