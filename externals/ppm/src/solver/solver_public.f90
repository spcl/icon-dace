!>
!! @file solver_public.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!! @author Florian Wilhelm <Florian.Wilhelm@kit.edu>
!
! Keywords: scales ppm solver parameter
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
!> provides naming parameters for the solver

MODULE solver_public

  IMPLICIT NONE

  !-------------------------------------------------------------------
  ! Parameters used in the way of enums
  !-------------------------------------------------------------------

  ! All solver types
  INTEGER, PARAMETER, PUBLIC :: CG_SOLVER        = 1
  INTEGER, PARAMETER, PUBLIC :: CHEBYSHEV_SOLVER = 2

  ! All preconditioner types
  INTEGER, PARAMETER, PUBLIC :: NONE_PRECOND     = 0
  INTEGER, PARAMETER, PUBLIC :: JACOBI_PRECOND   = 1
  INTEGER, PARAMETER, PUBLIC :: ILU0_PRECOND     = 2
  INTEGER, PARAMETER, PUBLIC :: SSOR_PRECOND     = 3
  INTEGER, PARAMETER, PUBLIC :: ICC_PRECOND      = 4
  INTEGER, PARAMETER, PUBLIC :: MICC_PRECOND     = 5

END MODULE solver_public
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
