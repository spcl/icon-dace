!>
!! @file solver_config_multi.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!
! Keywords: scales ppm solver configuration
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
#define filename 'solver_config_multi.f90'
!
! Check if stencil is defined
! Prototype is need to determine the right stencil, sp or dp
FUNCTION STENCIL_DEFINED(prototype) RESULT(ans)
  REAL(PREC), INTENT(IN) :: prototype
  LOGICAL :: ans

  ans = RANGE(prototype) == RANGE(STENCIL%shift) .AND. &
       ASSOCIATED(STENCIL%zonal) .AND. &
       ASSOCIATED(STENCIL%meridional) .AND. &
       ASSOCIATED(STENCIL%central) .AND. &
       ASSOCIATED(STENCIL%extent)

END FUNCTION STENCIL_DEFINED

! Set the matrix stencil for a symmetric 5-point-stencil
SUBROUTINE SET_STENCIL(zonal, meridional, central, ext)
  REAL(PREC), DIMENSION(:,:), TARGET, INTENT(IN) :: zonal, meridional, central
  TYPE(extent), DIMENSION(:), TARGET, INTENT(IN) :: ext

  STENCIL%zonal => zonal
  STENCIL%meridional => meridional
  STENCIL%central => central
  STENCIL%extent => ext

END SUBROUTINE SET_STENCIL

! This function applies the defined stencil or in other words
! the matrix vector multiplication
SUBROUTINE APPLY_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)
  INTEGER :: ib, il, jb, jl, i, j
  REAL(PREC), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  IF (.NOT. stencil_defined(0._precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         filename, __LINE__)
  ENDIF

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Define some aliases
  uf => STENCIL%zonal
  vf => STENCIL%meridional
  ff => STENCIL%central

  ! Apply matrix stencil
  DO j=jb,jl
    DO i=ib,il
      res_field(i,j) = ff(i,j)*field(i,j)                            &
           - uf(i,j)*field(i+1,j) - uf(i-1,j)*field(i-1,j) &
           - vf(i,j)*field(i,j+1) - vf(i,j-1)*field(i,j-1)
    ENDDO
  ENDDO

END SUBROUTINE APPLY_STENCIL

! This function applies the defined stencil or in other words
! the matrix vector multiplication but shifted with the defined value
SUBROUTINE APPLY_SHIFTED_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)
  INTEGER :: ib, il, jb, jl, i, j
  REAL(PREC), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  IF (.NOT. stencil_defined(0._precs)) THEN
    CALL abort_ppm("Stencil is not defined! Use set_stencil() first.", &
         filename, __LINE__)
  ENDIF

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Define some aliases
  uf => STENCIL%zonal
  vf => STENCIL%meridional
  ff => STENCIL%central

  ! Apply matrix stencil
  DO j=jb,jl
    DO i=ib,il
      res_field(i,j) = STENCIL%shift*field(i,j) - ff(i,j)*field(i,j) &
           + uf(i,j)*field(i+1,j) + uf(i-1,j)*field(i-1,j) &
           + vf(i,j)*field(i,j+1) + vf(i,j-1)*field(i,j-1)
    ENDDO
  ENDDO

END SUBROUTINE APPLY_SHIFTED_STENCIL

#undef filename
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
