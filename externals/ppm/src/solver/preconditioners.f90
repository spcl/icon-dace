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

!>
!! @file preconditioners_multi.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!
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
! Identity function, use this preconditioner for the non-precond. CG-method
SUBROUTINE identity_sp(r)
  REAL(sp), INTENT(INOUT) :: r(:,:)
  IF (SIZE(r, 1) > 0 .AND. SIZE(r, 2) > 0) THEN
    r(1, 1) = r(1, 1)
  END IF
END SUBROUTINE identity_sp

! Preparation for Jacobi preconditioner
SUBROUTINE prep_jacobi_sp(prototype)
  REAL(sp), INTENT(IN) :: prototype
  INTEGER :: ib, jb, il, jl, i, j

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ALLOCATE(jacobi_diag_sp(ib:il,jb:jl))

  DO j=jb,jl
    DO i=ib,il
      jacobi_diag_sp(i,j) = 1.0_precs/stencil_sp%central(i,j)
    ENDDO
  ENDDO

END SUBROUTINE prep_jacobi_sp

! Jacobi preconditioner
SUBROUTINE jacobi_sp(r)
  REAL(sp), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, jb, il, jl, i, j

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  DO j=jb,jl
    DO i=ib,il
      r(i,j) = r(i,j)*jacobi_diag_sp(i,j)
    ENDDO
  ENDDO

END SUBROUTINE jacobi_sp

! Preparation of the ILU(0) preconditioner
SUBROUTINE prep_ilu0_sp(prototype)
  REAL(sp), INTENT(in) :: prototype
  INTEGER :: ib, jb, il, jl, i, j
  REAL(sp), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ALLOCATE(ilu0_diag_sp(ib:il,jb:jl))

  ! Define some aliases
  uf => stencil_sp%zonal
  vf => stencil_sp%meridional
  ff => stencil_sp%central

  ! Calculate diagonal elements of ILU(0) factorization
  ilu0_diag_sp = 0.0_precs
  ilu0_diag_sp(ib,jb) = ff(ib,jb)
  DO i=ib+1,il
    ilu0_diag_sp(i,jb) = ff(i,jb) - uf(i-1,jb)**2/ilu0_diag_sp(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    ilu0_diag_sp(ib,j) = ff(ib,j) - vf(ib,j-1)**2/ilu0_diag_sp(ib,j-1)
  ENDDO
  DO j=jb+1,jl
    DO i=ib+1,il
      ilu0_diag_sp(i,j) = ff(i,j) - uf(i-1,j)**2/ilu0_diag_sp(i-1,j) &
           - vf(i,j-1)**2/ilu0_diag_sp(i,j-1)
    ENDDO
  ENDDO

END SUBROUTINE prep_ilu0_sp

! ILU(0) preconditioner
SUBROUTINE ilu0_sp(r)
  REAL(sp), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, jb, il, jl, i, j
  REAL(sp), DIMENSION(:,:), POINTER :: uf, vf

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Define some aliases
  uf => stencil_sp%zonal
  vf => stencil_sp%meridional

  ! The overall system to solve is As = r <-> (LU)s = r

  ! 1.) Solve system Ly = r
  DO i=ib+1,il
    r(i,jb) = r(i,jb) + uf(i-1,jb)/ilu0_diag_sp(i-1,jb)*r(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    r(ib,j) = r(ib,j) + vf(ib,j-1)/ilu0_diag_sp(ib,j-1)*r(ib,j-1)
  ENDDO
  DO j = jb+1,jl
    DO i = ib+1,il
      r(i,j) = r(i,j) + uf(i-1,j)/ilu0_diag_sp(i-1,j)*r(i-1,j) &
           + vf(i,j-1)/ilu0_diag_sp(i,j-1)*r(i,j-1)
    ENDDO
  ENDDO

  ! 2.) Solve system Us = y
  r(il,jl) = 1.0_precs/ilu0_diag_sp(il,jl)*r(il,jl)
  DO i=il-1,ib,-1
    r(i,jl) = 1.0_precs/ilu0_diag_sp(i,jl)*( r(i,jl) + uf(i,jl)*r(i+1,jl) )
  ENDDO
  DO j=jl-1,jb,-1
    r(il,j) = 1.0_precs/ilu0_diag_sp(il,j)*( r(il,j) + vf(il,j)*r(il,j+1) )
  ENDDO
  DO j = jl-1,jb,-1
    DO i = il-1,ib,-1
      r(i,j) = 1.0_precs/ilu0_diag_sp(i,j)*( r(i,j) + uf(i,j)*r(i+1,j)   &
           + vf(i,j)*r(i,j+1) )
    ENDDO
  ENDDO

END SUBROUTINE ilu0_sp

! Symmetric SOR preconditioner
SUBROUTINE SSOR_sp(r)

  REAL(sp), INTENT(INOUT) :: r(:,:)
  INTEGER :: i, j, ib, jb, il, jl
  REAL(sp):: ssorpar
  REAL(sp), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! Get SSOR_sp parameter from configuration
  ssorpar = REAL(get_ssor_param(), sp)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Define some aliases
  uf => stencil_sp%zonal
  vf => stencil_sp%meridional
  ff => stencil_sp%central

  !r(ib,jb) = r(ib,jb)
  DO i=ib+1,il
    r(i,jb) = r(i,jb) + ssorpar*uf(i-1,jb)/ff(i-1,jb)*r(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    r(ib,j) = r(ib,j) + ssorpar*vf(ib,j-1)/ff(ib,j-1)*r(ib,j-1)
  ENDDO
  DO j=jb+1,jl
    DO i=ib+1,il
      r(i,j) = r(i,j) + ssorpar*( uf(i-1,j)/ff(i-1,j)*r(i-1,j)  &
           + vf(i,j-1)/ff(i,j-1)*r(i,j-1) )
    ENDDO
  ENDDO

  r(il,jl) = r(il,jl)/ff(il,jl)
  DO i=il-1,ib,-1
    r(i,jl) = ( r(i,jl) + ssorpar*uf(i,jl)*r(i+1,jl) )/ff(i,jl)
  ENDDO
  DO j=jl-1,jb,-1
    r(il,j) = ( r(il,j) + ssorpar*vf(il,j)*r(il,j+1) )/ff(il,j)
  ENDDO
  DO j = jl-1,jb,-1
    DO i = il-1,ib,-1
      r(i,j) = ( r(i,j) + ssorpar*( uf(i,j)*r(i+1,j)           &
           + vf(i,j)*r(i,j+1) ) )/ff(i,j)
    ENDDO
  ENDDO

END SUBROUTINE SSOR_sp

! Calculates Matrix for Incomplete Cholesky Decomposition (icc_sp) with fillin p
! This routine makes use of a transformation from the barotropic stencil to a
! real matrix A with indices m_i, m_j. Variables s_i, s_j refer to the stencil
! coordinates, i, j are used according the actual context.
! The result L of L*L' ~= A is saved in ICC_C, ICC_W, ICC_S
SUBROUTINE prep_icc_sp(p,  prototype, modify_opt)
  INTEGER, INTENT(IN) :: p  ! FILL-IN
  REAL(sp), INTENT(IN) :: prototype
  LOGICAL, INTENT(IN), OPTIONAL :: modify_opt
  LOGICAL :: modify
  INTEGER :: ib, il, jb, jl, m_i, m_j, m_n, m_m, k
  REAL(sp) :: tmp, d, cs
  REAL(sp), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Define some aliases
  uf => stencil_sp%zonal
  vf => stencil_sp%meridional
  ff => stencil_sp%central

  ! Check and set optional arguments
  ! Check if modified icc_sp should be calculated
  IF ( PRESENT(modify_opt) ) THEN
    modify = modify_opt
  ELSE
    modify = .FALSE.
  ENDIF

  ! MAX(1,p-1) handles structurel exception for p=1
  ! ICC_C_sp = Center, _W = West, _S = South in stencil
  ALLOCATE(ICC_C_sp(ib:il,jb:jl), ICC_W_sp(MAX(1,p-1),ib:il,jb:jl), &
       ICC_S_sp(1+p,ib:il,jb:jl))

  ! Initialize
  ICC_C_sp = 1.0_precs
  ICC_W_sp = 0.0_precs
  ICC_S_sp = 0.0_precs

  m_n = extent_size(stencil_sp%extent)    ! dim of barotropic Matrix, m_ = Matrix
  m_m = extent_size(stencil_sp%extent(1)) ! m_m is dist between diag and outer diag

  DO m_i = 1,m_n

    ! start treating diagonal
    ! L(i,i) = sqrt( A(i,i) - sum(L(i,1:i-1).^2) ) in OCTAVE
    tmp = get_value_of_A(m_i,m_i)

    ! outer diagonal
    DO k = MAX(1,m_i-m_m),m_i-m_m+p
      tmp = tmp - get_value_of_L(m_i,k)**2
    ENDDO

    ! secondary diagonals
    DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
      tmp = tmp - get_value_of_L(m_i,k)**2
    ENDDO

    tmp = SQRT(tmp)
    CALL set_value_of_L(m_i, m_i, tmp)
    ! end treating diagonal

    ! start treating secondary diagonals
    DO m_j = m_i+1,m_i+MAX(1,p-1)
      IF ( m_j > m_n ) EXIT
      ! L(j,i) = ( A(j,i) - sum(L(j,1:i-1).*L(i,1:i-1)) )/L(i,i) in OCTAVE
      tmp = get_value_of_A(m_j, m_i)

      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      tmp = tmp / get_value_of_L(m_i,m_i)
      CALL set_value_of_L(m_j, m_i, tmp)
    ENDDO
    ! end treating secondary diagonals

    ! start treating outer diagonals
    DO m_j = m_i+m_m-p,m_i+m_m
      IF ( m_j > m_n ) EXIT
      ! L(j,i) = ( A(j,i) - sum(L(j,1:i-1).*L(i,1:i-1)) )/L(i,i) in OCTAVE
      tmp = get_value_of_A(m_j, m_i)

      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      tmp = tmp / get_value_of_L(m_i,m_i)
      CALL set_value_of_L(m_j, m_i, tmp)
    ENDDO
    ! end treating outer diagonals

    ! start calculate modified incomplete Cholesky
    IF (modify) THEN
      d = 0._precs ! holds defect in row to correct

      ! Left of diagonal
      ! start substract outer diagonal entry (m_i, m_i - m_m + p + 1)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - m_m + p + 1, k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - m_m + p + 1, k)
      ENDDO

      ! start substract secondary diagonal entry (m_i, m_i - p)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - MAX(2,p), k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - MAX(2,p), k)
      ENDDO

      ! Right of diagonal
      ! start substract outer diagonal entry (m_i, m_i + m_m - p - 1)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + m_m - p - 1, k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + m_m - p - 1, k)
      ENDDO

      ! start substract secondary diagonal entry (m_i, m_i + p)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + MAX(2,p), k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + MAX(2,p), k)
      ENDDO

      ! holds row sum with doubled diagonal entry
      cs = 2._precs * get_value_of_L(m_i,m_i)
      ! secondary diagonals
      DO m_j = m_i+1,m_i+MAX(1,p-1)
        cs = cs + get_value_of_L(m_j,m_i)
      ENDDO
      ! outer diagonals
      DO m_j = m_i+m_m-p,m_i+m_m
        cs = cs + get_value_of_L(m_j,m_i)
      ENDDO

      ! correct diagonal entry of L
      tmp = get_value_of_L(m_i, m_i)
      tmp = tmp + (-cs + SQRT(cs**2-4._precs*d)) / 2._precs
      CALL set_value_of_L(m_i, m_i, tmp)

    ENDIF

  ENDDO

  ! Calculate inverse of ICC_C_sp to avoid division afterwards
  ICC_C_sp = 1.0_precs/ICC_C_sp

CONTAINS

  ! Calculates stencil coordinates given a row in matrix A
  SUBROUTINE coord_by_row(row, i, j)
    INTEGER, INTENT(IN) :: row
    INTEGER, INTENT(OUT) :: i, j

    i = MOD(row - 1, m_m) + ib
    j = (row - 1) / m_m + jb

  END SUBROUTINE coord_by_row

  ! Get value in matrix L
  FUNCTION get_value_of_L(i,j) RESULT(res)
    INTEGER, INTENT(IN) :: i, j
    INTEGER :: d, s_i, s_j
    REAL(sp) :: res

    res = 0.0_precs

    ! return 0 if indices not in matrix (used if modify = True)
    IF ( i < 1 .OR. j < 1 .OR. i > m_n .OR. j > m_n ) RETURN

    ! diagonal
    IF ( i == j) THEN
      CALL coord_by_row(i, s_i, s_j)
      res = ICC_C_sp(s_i, s_j)
      RETURN
    ENDIF

    ! secondary diagonal
    d = i - j
    IF ( d >= 1 .AND. d <= MAX(1,p-1) ) THEN
      CALL coord_by_row(i, s_i, s_j)
      res = ICC_W_sp(d, s_i, s_j)
      RETURN
    ENDIF

    ! outer diagonal
    IF ( i - m_m > 0 ) THEN
      d = j - (i - m_m) + 1
      IF ( d >= 1 .AND. d <= p+1 ) THEN
        CALL coord_by_row(i, s_i, s_j)
        res = ICC_S_sp(d, s_i, s_j)
        RETURN
      ENDIF
    ENDIF

  END FUNCTION get_value_of_L

  ! Set value in matrix L
  SUBROUTINE set_value_of_L(i,j,v)
    INTEGER, INTENT(IN) :: i, j
    REAL(sp), INTENT(IN) :: v
    INTEGER :: d, s_i, s_j

    ! diagonal
    IF ( i == j) THEN
      CALL coord_by_row(i, s_i, s_j)
      ICC_C_sp(s_i, s_j) = v
      RETURN
    ENDIF

    ! secondary diagonal
    d = i - j
    IF ( d >= 1 .AND. d <= MAX(1,p-1) ) THEN
      CALL coord_by_row(i, s_i, s_j)
      ICC_W_sp(d, s_i, s_j) = v
      RETURN
    ENDIF

    ! outer diagonal
    IF ( i - m_m > 0 ) THEN
      d = j - (i - m_m) + 1
      IF ( d >= 1 .AND. d <= p+1 ) THEN
        CALL coord_by_row(i, s_i, s_j)
        ICC_S_sp(d, s_i, s_j) = v
        RETURN
      ENDIF
    ENDIF

  END SUBROUTINE set_value_of_L

  ! Get value of virtual barotropic matrix A
  FUNCTION get_value_of_A(i,j) RESULT(res)
    INTEGER, INTENT(IN) :: i, j
    INTEGER :: s_i, s_j
    REAL(sp) :: res

    res = 0.0_precs

    IF ( i == j) THEN ! diagonal
      CALL coord_by_row(i, s_i, s_j)
      res = ff(s_i, s_j)
    ELSEIF ( j-i == 1 ) THEN ! right secondary diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_i >= il) RETURN ! hit boundary
      res = -uf(s_i,s_j)
    ELSEIF ( j-i == -1 ) THEN ! left secondary diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_i <= ib) RETURN ! hit boundary
      res = -uf(s_i-1,s_j)
    ELSEIF ( j-i == m_m ) THEN ! right outer diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_j >= jl) RETURN ! hit boundary
      res = -vf(s_i,s_j)
    ELSEIF ( j-i == -m_m ) THEN ! left outer diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_j <= jb) RETURN ! hit boundary
      res = -vf(s_i,s_j-1)
    ENDIF

  END FUNCTION get_value_of_A

END SUBROUTINE prep_icc_sp

! icc_sp(p) preconditioner
! Attention: This function silently assumes that r is an ordered real
! (not NaN nor +-Inf) outside of r(ib:il,jb:jl)
SUBROUTINE icc_sp(r)
  REAL(sp), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, il, jb, jl, i, j, k, p, s_l, w_l

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Set halos to zero would be more secure at this point but affects
  ! the performance negatively. This should be done outside of this.
  !CALL clear_halos(r, stencil%extent)

  ! get the fill-in p of icc_sp
  p = SIZE(ICC_S_sp,1) - 1
  s_l = p + 1 ! length of south diagonal
  w_l = MAX(1,p-1) ! length of west diagonal

  ! The overall system to solve is As = r <-> (LL')s = r

  ! 1.) Solve system Ly = r

  ! first column
  r(ib,jb) = r(ib,jb)*ICC_C_sp(ib,jb)
  DO i=ib+1,il
    r(i,jb) = ( r(i,jb) - ICC_W_sp(1,i,jb)*r(i-1,jb) )*ICC_C_sp(i,jb)
  ENDDO

  ! all other columns
  DO j = jb+1,jl
    DO i = ib,il
      DO k = 1,s_l
        r(i,j) = r(i,j) - ICC_S_sp(k,i,j)*r(i+k-1,j-1)
      ENDDO
      DO k = 1,w_l
        r(i,j) = r(i,j) - ICC_W_sp(k,i,j)*r(i-k,j)
      ENDDO
      r(i,j) = r(i,j)*ICC_C_sp(i,j)
    ENDDO
  ENDDO

  ! 2.) Solve system L's = y
  ! Backward substitution is done on a column level instead of row-wise

  ! all but first columns
  DO j = jl,jb+1,-1
    DO i = il,ib,-1
      r(i,j) = r(i,j)*ICC_C_sp(i,j)
      DO k = w_l,1,-1
        r(i-k,j) = r(i-k,j) - r(i,j)*ICC_W_sp(k,i,j) ! saxpy
      ENDDO
      DO k = s_l,1,-1
        r(i+k-1,j-1) = r(i+k-1,j-1) - r(i,j)*ICC_S_sp(k,i,j) ! saxpy
      ENDDO
    ENDDO
  ENDDO

  ! first column
  DO i = il,ib+1,-1
    r(i,jb) = r(i,jb)*ICC_C_sp(i,jb)
    r(i-1,jb) = r(i-1,jb) - r(i,jb)*ICC_W_sp(1,i,jb) ! saxpy
  ENDDO
  r(ib,jb) = r(ib,jb)*ICC_C_sp(ib,jb)

END SUBROUTINE icc_sp

! Prepare MICC(p) preconditioner
! Modified Incomplete Cholesky with fill-in p. The modification is that the
! diagonal entries of icc_sp(p) get modified in the way that the row sums of
! L*L' equal the row sums of A.
! MICC(p) only differs from icc_sp(p) in the preparation step. The application
! of MICC(p) is done by calling the subroutine icc.
SUBROUTINE prep_micc_sp(p, prototype)
  INTEGER, INTENT(IN) :: p  ! FILL-IN
  REAL(sp), INTENT(IN) :: prototype
  CALL prep_icc(p, prototype, .TRUE.)
END SUBROUTINE prep_micc_sp

! Check if a certain preconditioner is already prepared
FUNCTION precond_prepared_sp(preconditioner, prototype) RESULT(ans)
  INTEGER, INTENT(IN) :: preconditioner
  REAL(sp), INTENT(IN) :: prototype
  LOGICAL :: ans
  CHARACTER(len=25), PARAMETER :: msg_prefix = 'Selected preconditioner ['
  CHARACTER(len=17), PARAMETER :: msg_suffix = '] does not exist!'
  CHARACTER(len=LEN(msg_prefix) + LEN(msg_suffix) + 11) :: msg

  ans = KIND(prototype) == sp
  SELECT CASE (preconditioner)
  CASE (NONE_PRECOND)
    ans = .TRUE.
  CASE (JACOBI_PRECOND)
    ans = ALLOCATED(jacobi_diag_sp)
  CASE (ILU0_PRECOND)
    ans = ALLOCATED(ilu0_diag_sp)
  CASE (SSOR_PRECOND)
    ans = .NOT. ieee_is_nan(config%ssor_param)
  CASE (ICC_PRECOND,MICC_PRECOND)
    ans = ALLOCATED(ICC_C_sp) .AND. ALLOCATED(ICC_W_sp) .AND. ALLOCATED(ICC_S_sp)
  CASE DEFAULT
    WRITE (msg, '(a,i0,a)') msg_prefix, preconditioner, msg_suffix
    CALL abort_ppm(msg, 'preconditioners_multi.f90', 651)
  END SELECT

END FUNCTION precond_prepared_sp

! Jacobi preconditioned matrix stencil operation
SUBROUTINE jacobi_precond_stencil_sp(field, res_field)

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)

  ! Apply stencil first
  IF (.NOT. precond_prepared(JACOBI_PRECOND, 0._precs)) THEN
    CALL prep_jacobi(0._precs)
  ENDIF
  CALL apply_stencil(field, res_field)

  ! Apply Jacobi preconditioner
  CALL jacobi(res_field)

END SUBROUTINE jacobi_precond_stencil_sp

! Jacobi preconditioned matrix stencil operation shifted by a value
SUBROUTINE jacobi_precond_shifted_stencil_sp(field, res_field)

  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Apply stencil first
  IF (.NOT. precond_prepared(JACOBI_PRECOND, 0._precs)) THEN
    CALL prep_jacobi(0._precs)
  ENDIF
  CALL apply_stencil(field, res_field)

  ! Apply Jacobi preconditioner
  CALL jacobi(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_sp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE jacobi_precond_shifted_stencil_sp

! ILU(0) preconditioned matrix stencil operation
SUBROUTINE ilu0_precond_stencil_sp(field, res_field)
  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ilu0_sp preconditioner
  IF (.NOT. precond_prepared(ILU0_PRECOND, 0._precs)) THEN
    CALL prep_ilu0(0._precs)
  ENDIF
  CALL ilu0(res_field)

END SUBROUTINE ilu0_precond_stencil_sp

! ILU(0) preconditioned matrix stencil operation shifted by a value
SUBROUTINE ilu0_precond_shifted_stencil_sp(field, res_field)
  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ilu0_sp preconditioner
  IF (.NOT. precond_prepared(ILU0_PRECOND, 0._precs)) THEN
    CALL prep_ilu0(0._precs)
  ENDIF
  CALL ilu0(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_sp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE ilu0_precond_shifted_stencil_sp

! SSOR_sp preconditioned matrix stencil operation
SUBROUTINE ssor_precond_stencil_sp(field, res_field)
  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply SSOR_sp preconditioner
  IF (.NOT. precond_prepared(SSOR_PRECOND, 0._precs)) THEN
    CALL abort_ppm("No SSOR parameter provided!", 'preconditioners_multi.f90', 775)
  ENDIF
  CALL ssor(res_field)

END SUBROUTINE ssor_precond_stencil_sp

! SSOR_sp preconditioned matrix stencil operation shifted by a value
SUBROUTINE ssor_precond_shifted_stencil_sp(field, res_field)
  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply SSOR_sp preconditioner
  IF (.NOT. precond_prepared(SSOR_PRECOND, 0._precs)) THEN
    CALL abort_ppm("No SSOR parameter provided!", 'preconditioners_multi.f90', 801)
  ENDIF
  CALL ssor(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_sp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE ssor_precond_shifted_stencil_sp

! icc_sp(p) preconditioned matrix stencil operation
SUBROUTINE icc_precond_stencil_sp(field, res_field)
  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply icc_sp preconditioner
  IF (.NOT. precond_prepared(ICC_PRECOND, 0._precs)) THEN
    CALL prep_icc(config%icc_param, 0._precs)
  ENDIF
  CALL icc(res_field)

END SUBROUTINE icc_precond_stencil_sp

! icc_sp(p) preconditioned matrix stencil operation shifted by a value
SUBROUTINE icc_precond_shifted_stencil_sp(field, res_field)
  REAL(sp), INTENT(IN) :: field(:,:)
  REAL(sp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! sp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = sp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_sp%extent(1))
  jb = extent_start(stencil_sp%extent(2))
  il = extent_end(stencil_sp%extent(1))
  jl = extent_end(stencil_sp%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply icc_sp preconditioner
  IF (.NOT. precond_prepared(ICC_PRECOND, 0._precs)) THEN
    CALL prep_icc(config%icc_param, 0._precs)
  ENDIF
  CALL icc(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_sp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE icc_precond_shifted_stencil_sp
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
!>
!! @file preconditioners_multi.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!
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
! Identity function, use this preconditioner for the non-precond. CG-method
SUBROUTINE identity_dp(r)
  REAL(dp), INTENT(INOUT) :: r(:,:)
  IF (SIZE(r, 1) > 0 .AND. SIZE(r, 2) > 0) THEN
    r(1, 1) = r(1, 1)
  END IF
END SUBROUTINE identity_dp

! Preparation for Jacobi preconditioner
SUBROUTINE prep_jacobi_dp(prototype)
  REAL(dp), INTENT(IN) :: prototype
  INTEGER :: ib, jb, il, jl, i, j

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ALLOCATE(jacobi_diag_dp(ib:il,jb:jl))

  DO j=jb,jl
    DO i=ib,il
      jacobi_diag_dp(i,j) = 1.0_precs/stencil_dp%central(i,j)
    ENDDO
  ENDDO

END SUBROUTINE prep_jacobi_dp

! Jacobi preconditioner
SUBROUTINE jacobi_dp(r)
  REAL(dp), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, jb, il, jl, i, j

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  DO j=jb,jl
    DO i=ib,il
      r(i,j) = r(i,j)*jacobi_diag_dp(i,j)
    ENDDO
  ENDDO

END SUBROUTINE jacobi_dp

! Preparation of the ILU(0) preconditioner
SUBROUTINE prep_ilu0_dp(prototype)
  REAL(dp), INTENT(in) :: prototype
  INTEGER :: ib, jb, il, jl, i, j
  REAL(dp), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ALLOCATE(ilu0_diag_dp(ib:il,jb:jl))

  ! Define some aliases
  uf => stencil_dp%zonal
  vf => stencil_dp%meridional
  ff => stencil_dp%central

  ! Calculate diagonal elements of ILU(0) factorization
  ilu0_diag_dp = 0.0_precs
  ilu0_diag_dp(ib,jb) = ff(ib,jb)
  DO i=ib+1,il
    ilu0_diag_dp(i,jb) = ff(i,jb) - uf(i-1,jb)**2/ilu0_diag_dp(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    ilu0_diag_dp(ib,j) = ff(ib,j) - vf(ib,j-1)**2/ilu0_diag_dp(ib,j-1)
  ENDDO
  DO j=jb+1,jl
    DO i=ib+1,il
      ilu0_diag_dp(i,j) = ff(i,j) - uf(i-1,j)**2/ilu0_diag_dp(i-1,j) &
           - vf(i,j-1)**2/ilu0_diag_dp(i,j-1)
    ENDDO
  ENDDO

END SUBROUTINE prep_ilu0_dp

! ILU(0) preconditioner
SUBROUTINE ilu0_dp(r)
  REAL(dp), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, jb, il, jl, i, j
  REAL(dp), DIMENSION(:,:), POINTER :: uf, vf

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Define some aliases
  uf => stencil_dp%zonal
  vf => stencil_dp%meridional

  ! The overall system to solve is As = r <-> (LU)s = r

  ! 1.) Solve system Ly = r
  DO i=ib+1,il
    r(i,jb) = r(i,jb) + uf(i-1,jb)/ilu0_diag_dp(i-1,jb)*r(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    r(ib,j) = r(ib,j) + vf(ib,j-1)/ilu0_diag_dp(ib,j-1)*r(ib,j-1)
  ENDDO
  DO j = jb+1,jl
    DO i = ib+1,il
      r(i,j) = r(i,j) + uf(i-1,j)/ilu0_diag_dp(i-1,j)*r(i-1,j) &
           + vf(i,j-1)/ilu0_diag_dp(i,j-1)*r(i,j-1)
    ENDDO
  ENDDO

  ! 2.) Solve system Us = y
  r(il,jl) = 1.0_precs/ilu0_diag_dp(il,jl)*r(il,jl)
  DO i=il-1,ib,-1
    r(i,jl) = 1.0_precs/ilu0_diag_dp(i,jl)*( r(i,jl) + uf(i,jl)*r(i+1,jl) )
  ENDDO
  DO j=jl-1,jb,-1
    r(il,j) = 1.0_precs/ilu0_diag_dp(il,j)*( r(il,j) + vf(il,j)*r(il,j+1) )
  ENDDO
  DO j = jl-1,jb,-1
    DO i = il-1,ib,-1
      r(i,j) = 1.0_precs/ilu0_diag_dp(i,j)*( r(i,j) + uf(i,j)*r(i+1,j)   &
           + vf(i,j)*r(i,j+1) )
    ENDDO
  ENDDO

END SUBROUTINE ilu0_dp

! Symmetric SOR preconditioner
SUBROUTINE SSOR_dp(r)

  REAL(dp), INTENT(INOUT) :: r(:,:)
  INTEGER :: i, j, ib, jb, il, jl
  REAL(dp):: ssorpar
  REAL(dp), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! Get SSOR_dp parameter from configuration
  ssorpar = REAL(get_ssor_param(), dp)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Define some aliases
  uf => stencil_dp%zonal
  vf => stencil_dp%meridional
  ff => stencil_dp%central

  !r(ib,jb) = r(ib,jb)
  DO i=ib+1,il
    r(i,jb) = r(i,jb) + ssorpar*uf(i-1,jb)/ff(i-1,jb)*r(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    r(ib,j) = r(ib,j) + ssorpar*vf(ib,j-1)/ff(ib,j-1)*r(ib,j-1)
  ENDDO
  DO j=jb+1,jl
    DO i=ib+1,il
      r(i,j) = r(i,j) + ssorpar*( uf(i-1,j)/ff(i-1,j)*r(i-1,j)  &
           + vf(i,j-1)/ff(i,j-1)*r(i,j-1) )
    ENDDO
  ENDDO

  r(il,jl) = r(il,jl)/ff(il,jl)
  DO i=il-1,ib,-1
    r(i,jl) = ( r(i,jl) + ssorpar*uf(i,jl)*r(i+1,jl) )/ff(i,jl)
  ENDDO
  DO j=jl-1,jb,-1
    r(il,j) = ( r(il,j) + ssorpar*vf(il,j)*r(il,j+1) )/ff(il,j)
  ENDDO
  DO j = jl-1,jb,-1
    DO i = il-1,ib,-1
      r(i,j) = ( r(i,j) + ssorpar*( uf(i,j)*r(i+1,j)           &
           + vf(i,j)*r(i,j+1) ) )/ff(i,j)
    ENDDO
  ENDDO

END SUBROUTINE SSOR_dp

! Calculates Matrix for Incomplete Cholesky Decomposition (icc_dp) with fillin p
! This routine makes use of a transformation from the barotropic stencil to a
! real matrix A with indices m_i, m_j. Variables s_i, s_j refer to the stencil
! coordinates, i, j are used according the actual context.
! The result L of L*L' ~= A is saved in ICC_C, ICC_W, ICC_S
SUBROUTINE prep_icc_dp(p,  prototype, modify_opt)
  INTEGER, INTENT(IN) :: p  ! FILL-IN
  REAL(dp), INTENT(IN) :: prototype
  LOGICAL, INTENT(IN), OPTIONAL :: modify_opt
  LOGICAL :: modify
  INTEGER :: ib, il, jb, jl, m_i, m_j, m_n, m_m, k
  REAL(dp) :: tmp, d, cs
  REAL(dp), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Define some aliases
  uf => stencil_dp%zonal
  vf => stencil_dp%meridional
  ff => stencil_dp%central

  ! Check and set optional arguments
  ! Check if modified icc_dp should be calculated
  IF ( PRESENT(modify_opt) ) THEN
    modify = modify_opt
  ELSE
    modify = .FALSE.
  ENDIF

  ! MAX(1,p-1) handles structurel exception for p=1
  ! ICC_C_dp = Center, _W = West, _S = South in stencil
  ALLOCATE(ICC_C_dp(ib:il,jb:jl), ICC_W_dp(MAX(1,p-1),ib:il,jb:jl), &
       ICC_S_dp(1+p,ib:il,jb:jl))

  ! Initialize
  ICC_C_dp = 1.0_precs
  ICC_W_dp = 0.0_precs
  ICC_S_dp = 0.0_precs

  m_n = extent_size(stencil_dp%extent)    ! dim of barotropic Matrix, m_ = Matrix
  m_m = extent_size(stencil_dp%extent(1)) ! m_m is dist between diag and outer diag

  DO m_i = 1,m_n

    ! start treating diagonal
    ! L(i,i) = sqrt( A(i,i) - sum(L(i,1:i-1).^2) ) in OCTAVE
    tmp = get_value_of_A(m_i,m_i)

    ! outer diagonal
    DO k = MAX(1,m_i-m_m),m_i-m_m+p
      tmp = tmp - get_value_of_L(m_i,k)**2
    ENDDO

    ! secondary diagonals
    DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
      tmp = tmp - get_value_of_L(m_i,k)**2
    ENDDO

    tmp = SQRT(tmp)
    CALL set_value_of_L(m_i, m_i, tmp)
    ! end treating diagonal

    ! start treating secondary diagonals
    DO m_j = m_i+1,m_i+MAX(1,p-1)
      IF ( m_j > m_n ) EXIT
      ! L(j,i) = ( A(j,i) - sum(L(j,1:i-1).*L(i,1:i-1)) )/L(i,i) in OCTAVE
      tmp = get_value_of_A(m_j, m_i)

      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      tmp = tmp / get_value_of_L(m_i,m_i)
      CALL set_value_of_L(m_j, m_i, tmp)
    ENDDO
    ! end treating secondary diagonals

    ! start treating outer diagonals
    DO m_j = m_i+m_m-p,m_i+m_m
      IF ( m_j > m_n ) EXIT
      ! L(j,i) = ( A(j,i) - sum(L(j,1:i-1).*L(i,1:i-1)) )/L(i,i) in OCTAVE
      tmp = get_value_of_A(m_j, m_i)

      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        tmp = tmp - get_value_of_L(m_j,k) * get_value_of_L(m_i,k)
      ENDDO

      tmp = tmp / get_value_of_L(m_i,m_i)
      CALL set_value_of_L(m_j, m_i, tmp)
    ENDDO
    ! end treating outer diagonals

    ! start calculate modified incomplete Cholesky
    IF (modify) THEN
      d = 0._precs ! holds defect in row to correct

      ! Left of diagonal
      ! start substract outer diagonal entry (m_i, m_i - m_m + p + 1)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - m_m + p + 1, k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - m_m + p + 1, k)
      ENDDO

      ! start substract secondary diagonal entry (m_i, m_i - p)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - MAX(2,p), k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i - MAX(2,p), k)
      ENDDO

      ! Right of diagonal
      ! start substract outer diagonal entry (m_i, m_i + m_m - p - 1)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + m_m - p - 1, k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + m_m - p - 1, k)
      ENDDO

      ! start substract secondary diagonal entry (m_i, m_i + p)
      ! outer diagonal
      DO k = MAX(1,m_i-m_m),m_i-m_m+p
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + MAX(2,p), k)
      ENDDO

      ! secondary diagonals
      DO k = MAX(1,m_i-MAX(1,p-1)),m_i-1
        d = d + get_value_of_L(m_i,k) &
             * get_value_of_L(m_i + MAX(2,p), k)
      ENDDO

      ! holds row sum with doubled diagonal entry
      cs = 2._precs * get_value_of_L(m_i,m_i)
      ! secondary diagonals
      DO m_j = m_i+1,m_i+MAX(1,p-1)
        cs = cs + get_value_of_L(m_j,m_i)
      ENDDO
      ! outer diagonals
      DO m_j = m_i+m_m-p,m_i+m_m
        cs = cs + get_value_of_L(m_j,m_i)
      ENDDO

      ! correct diagonal entry of L
      tmp = get_value_of_L(m_i, m_i)
      tmp = tmp + (-cs + SQRT(cs**2-4._precs*d)) / 2._precs
      CALL set_value_of_L(m_i, m_i, tmp)

    ENDIF

  ENDDO

  ! Calculate inverse of ICC_C_dp to avoid division afterwards
  ICC_C_dp = 1.0_precs/ICC_C_dp

CONTAINS

  ! Calculates stencil coordinates given a row in matrix A
  SUBROUTINE coord_by_row(row, i, j)
    INTEGER, INTENT(IN) :: row
    INTEGER, INTENT(OUT) :: i, j

    i = MOD(row - 1, m_m) + ib
    j = (row - 1) / m_m + jb

  END SUBROUTINE coord_by_row

  ! Get value in matrix L
  FUNCTION get_value_of_L(i,j) RESULT(res)
    INTEGER, INTENT(IN) :: i, j
    INTEGER :: d, s_i, s_j
    REAL(dp) :: res

    res = 0.0_precs

    ! return 0 if indices not in matrix (used if modify = True)
    IF ( i < 1 .OR. j < 1 .OR. i > m_n .OR. j > m_n ) RETURN

    ! diagonal
    IF ( i == j) THEN
      CALL coord_by_row(i, s_i, s_j)
      res = ICC_C_dp(s_i, s_j)
      RETURN
    ENDIF

    ! secondary diagonal
    d = i - j
    IF ( d >= 1 .AND. d <= MAX(1,p-1) ) THEN
      CALL coord_by_row(i, s_i, s_j)
      res = ICC_W_dp(d, s_i, s_j)
      RETURN
    ENDIF

    ! outer diagonal
    IF ( i - m_m > 0 ) THEN
      d = j - (i - m_m) + 1
      IF ( d >= 1 .AND. d <= p+1 ) THEN
        CALL coord_by_row(i, s_i, s_j)
        res = ICC_S_dp(d, s_i, s_j)
        RETURN
      ENDIF
    ENDIF

  END FUNCTION get_value_of_L

  ! Set value in matrix L
  SUBROUTINE set_value_of_L(i,j,v)
    INTEGER, INTENT(IN) :: i, j
    REAL(dp), INTENT(IN) :: v
    INTEGER :: d, s_i, s_j

    ! diagonal
    IF ( i == j) THEN
      CALL coord_by_row(i, s_i, s_j)
      ICC_C_dp(s_i, s_j) = v
      RETURN
    ENDIF

    ! secondary diagonal
    d = i - j
    IF ( d >= 1 .AND. d <= MAX(1,p-1) ) THEN
      CALL coord_by_row(i, s_i, s_j)
      ICC_W_dp(d, s_i, s_j) = v
      RETURN
    ENDIF

    ! outer diagonal
    IF ( i - m_m > 0 ) THEN
      d = j - (i - m_m) + 1
      IF ( d >= 1 .AND. d <= p+1 ) THEN
        CALL coord_by_row(i, s_i, s_j)
        ICC_S_dp(d, s_i, s_j) = v
        RETURN
      ENDIF
    ENDIF

  END SUBROUTINE set_value_of_L

  ! Get value of virtual barotropic matrix A
  FUNCTION get_value_of_A(i,j) RESULT(res)
    INTEGER, INTENT(IN) :: i, j
    INTEGER :: s_i, s_j
    REAL(dp) :: res

    res = 0.0_precs

    IF ( i == j) THEN ! diagonal
      CALL coord_by_row(i, s_i, s_j)
      res = ff(s_i, s_j)
    ELSEIF ( j-i == 1 ) THEN ! right secondary diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_i >= il) RETURN ! hit boundary
      res = -uf(s_i,s_j)
    ELSEIF ( j-i == -1 ) THEN ! left secondary diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_i <= ib) RETURN ! hit boundary
      res = -uf(s_i-1,s_j)
    ELSEIF ( j-i == m_m ) THEN ! right outer diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_j >= jl) RETURN ! hit boundary
      res = -vf(s_i,s_j)
    ELSEIF ( j-i == -m_m ) THEN ! left outer diagonal
      CALL coord_by_row(i, s_i, s_j)
      IF (s_j <= jb) RETURN ! hit boundary
      res = -vf(s_i,s_j-1)
    ENDIF

  END FUNCTION get_value_of_A

END SUBROUTINE prep_icc_dp

! icc_dp(p) preconditioner
! Attention: This function silently assumes that r is an ordered real
! (not NaN nor +-Inf) outside of r(ib:il,jb:jl)
SUBROUTINE icc_dp(r)
  REAL(dp), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, il, jb, jl, i, j, k, p, s_l, w_l

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Set halos to zero would be more secure at this point but affects
  ! the performance negatively. This should be done outside of this.
  !CALL clear_halos(r, stencil%extent)

  ! get the fill-in p of icc_dp
  p = SIZE(ICC_S_dp,1) - 1
  s_l = p + 1 ! length of south diagonal
  w_l = MAX(1,p-1) ! length of west diagonal

  ! The overall system to solve is As = r <-> (LL')s = r

  ! 1.) Solve system Ly = r

  ! first column
  r(ib,jb) = r(ib,jb)*ICC_C_dp(ib,jb)
  DO i=ib+1,il
    r(i,jb) = ( r(i,jb) - ICC_W_dp(1,i,jb)*r(i-1,jb) )*ICC_C_dp(i,jb)
  ENDDO

  ! all other columns
  DO j = jb+1,jl
    DO i = ib,il
      DO k = 1,s_l
        r(i,j) = r(i,j) - ICC_S_dp(k,i,j)*r(i+k-1,j-1)
      ENDDO
      DO k = 1,w_l
        r(i,j) = r(i,j) - ICC_W_dp(k,i,j)*r(i-k,j)
      ENDDO
      r(i,j) = r(i,j)*ICC_C_dp(i,j)
    ENDDO
  ENDDO

  ! 2.) Solve system L's = y
  ! Backward substitution is done on a column level instead of row-wise

  ! all but first columns
  DO j = jl,jb+1,-1
    DO i = il,ib,-1
      r(i,j) = r(i,j)*ICC_C_dp(i,j)
      DO k = w_l,1,-1
        r(i-k,j) = r(i-k,j) - r(i,j)*ICC_W_dp(k,i,j) ! saxpy
      ENDDO
      DO k = s_l,1,-1
        r(i+k-1,j-1) = r(i+k-1,j-1) - r(i,j)*ICC_S_dp(k,i,j) ! saxpy
      ENDDO
    ENDDO
  ENDDO

  ! first column
  DO i = il,ib+1,-1
    r(i,jb) = r(i,jb)*ICC_C_dp(i,jb)
    r(i-1,jb) = r(i-1,jb) - r(i,jb)*ICC_W_dp(1,i,jb) ! saxpy
  ENDDO
  r(ib,jb) = r(ib,jb)*ICC_C_dp(ib,jb)

END SUBROUTINE icc_dp

! Prepare MICC(p) preconditioner
! Modified Incomplete Cholesky with fill-in p. The modification is that the
! diagonal entries of icc_dp(p) get modified in the way that the row sums of
! L*L' equal the row sums of A.
! MICC(p) only differs from icc_dp(p) in the preparation step. The application
! of MICC(p) is done by calling the subroutine icc.
SUBROUTINE prep_micc_dp(p, prototype)
  INTEGER, INTENT(IN) :: p  ! FILL-IN
  REAL(dp), INTENT(IN) :: prototype
  CALL prep_icc(p, prototype, .TRUE.)
END SUBROUTINE prep_micc_dp

! Check if a certain preconditioner is already prepared
FUNCTION precond_prepared_dp(preconditioner, prototype) RESULT(ans)
  INTEGER, INTENT(IN) :: preconditioner
  REAL(dp), INTENT(IN) :: prototype
  LOGICAL :: ans
  CHARACTER(len=25), PARAMETER :: msg_prefix = 'Selected preconditioner ['
  CHARACTER(len=17), PARAMETER :: msg_suffix = '] does not exist!'
  CHARACTER(len=LEN(msg_prefix) + LEN(msg_suffix) + 11) :: msg

  ans = KIND(prototype) == dp
  SELECT CASE (preconditioner)
  CASE (NONE_PRECOND)
    ans = .TRUE.
  CASE (JACOBI_PRECOND)
    ans = ALLOCATED(jacobi_diag_dp)
  CASE (ILU0_PRECOND)
    ans = ALLOCATED(ilu0_diag_dp)
  CASE (SSOR_PRECOND)
    ans = .NOT. ieee_is_nan(config%ssor_param)
  CASE (ICC_PRECOND,MICC_PRECOND)
    ans = ALLOCATED(ICC_C_dp) .AND. ALLOCATED(ICC_W_dp) .AND. ALLOCATED(ICC_S_dp)
  CASE DEFAULT
    WRITE (msg, '(a,i0,a)') msg_prefix, preconditioner, msg_suffix
    CALL abort_ppm(msg, 'preconditioners_multi.f90', 651)
  END SELECT

END FUNCTION precond_prepared_dp

! Jacobi preconditioned matrix stencil operation
SUBROUTINE jacobi_precond_stencil_dp(field, res_field)

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)

  ! Apply stencil first
  IF (.NOT. precond_prepared(JACOBI_PRECOND, 0._precs)) THEN
    CALL prep_jacobi(0._precs)
  ENDIF
  CALL apply_stencil(field, res_field)

  ! Apply Jacobi preconditioner
  CALL jacobi(res_field)

END SUBROUTINE jacobi_precond_stencil_dp

! Jacobi preconditioned matrix stencil operation shifted by a value
SUBROUTINE jacobi_precond_shifted_stencil_dp(field, res_field)

  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Apply stencil first
  IF (.NOT. precond_prepared(JACOBI_PRECOND, 0._precs)) THEN
    CALL prep_jacobi(0._precs)
  ENDIF
  CALL apply_stencil(field, res_field)

  ! Apply Jacobi preconditioner
  CALL jacobi(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_dp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE jacobi_precond_shifted_stencil_dp

! ILU(0) preconditioned matrix stencil operation
SUBROUTINE ilu0_precond_stencil_dp(field, res_field)
  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ilu0_dp preconditioner
  IF (.NOT. precond_prepared(ILU0_PRECOND, 0._precs)) THEN
    CALL prep_ilu0(0._precs)
  ENDIF
  CALL ilu0(res_field)

END SUBROUTINE ilu0_precond_stencil_dp

! ILU(0) preconditioned matrix stencil operation shifted by a value
SUBROUTINE ilu0_precond_shifted_stencil_dp(field, res_field)
  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ilu0_dp preconditioner
  IF (.NOT. precond_prepared(ILU0_PRECOND, 0._precs)) THEN
    CALL prep_ilu0(0._precs)
  ENDIF
  CALL ilu0(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_dp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE ilu0_precond_shifted_stencil_dp

! SSOR_dp preconditioned matrix stencil operation
SUBROUTINE ssor_precond_stencil_dp(field, res_field)
  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply SSOR_dp preconditioner
  IF (.NOT. precond_prepared(SSOR_PRECOND, 0._precs)) THEN
    CALL abort_ppm("No SSOR parameter provided!", 'preconditioners_multi.f90', 775)
  ENDIF
  CALL ssor(res_field)

END SUBROUTINE ssor_precond_stencil_dp

! SSOR_dp preconditioned matrix stencil operation shifted by a value
SUBROUTINE ssor_precond_shifted_stencil_dp(field, res_field)
  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply SSOR_dp preconditioner
  IF (.NOT. precond_prepared(SSOR_PRECOND, 0._precs)) THEN
    CALL abort_ppm("No SSOR parameter provided!", 'preconditioners_multi.f90', 801)
  ENDIF
  CALL ssor(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_dp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE ssor_precond_shifted_stencil_dp

! icc_dp(p) preconditioned matrix stencil operation
SUBROUTINE icc_precond_stencil_dp(field, res_field)
  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply icc_dp preconditioner
  IF (.NOT. precond_prepared(ICC_PRECOND, 0._precs)) THEN
    CALL prep_icc(config%icc_param, 0._precs)
  ENDIF
  CALL icc(res_field)

END SUBROUTINE icc_precond_stencil_dp

! icc_dp(p) preconditioned matrix stencil operation shifted by a value
SUBROUTINE icc_precond_shifted_stencil_dp(field, res_field)
  REAL(dp), INTENT(IN) :: field(:,:)
  REAL(dp), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! dp suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = dp

  ! Define some variables for the ranges of the fields
  ib = extent_start(stencil_dp%extent(1))
  jb = extent_start(stencil_dp%extent(2))
  il = extent_end(stencil_dp%extent(1))
  jl = extent_end(stencil_dp%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply icc_dp preconditioner
  IF (.NOT. precond_prepared(ICC_PRECOND, 0._precs)) THEN
    CALL prep_icc(config%icc_param, 0._precs)
  ENDIF
  CALL icc(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = stencil_dp%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE icc_precond_shifted_stencil_dp
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:

END MODULE preconditioners
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
