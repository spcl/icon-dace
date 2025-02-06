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
#define filename 'preconditioners_multi.f90'
!
! Identity function, use this preconditioner for the non-precond. CG-method
SUBROUTINE IDENTITY(r)
  REAL(PREC), INTENT(INOUT) :: r(:,:)
  IF (SIZE(r, 1) > 0 .AND. SIZE(r, 2) > 0) THEN
    r(1, 1) = r(1, 1)
  END IF
END SUBROUTINE IDENTITY

! Preparation for Jacobi preconditioner
SUBROUTINE PREP_JACOBI(prototype)
  REAL(PREC), INTENT(IN) :: prototype
  INTEGER :: ib, jb, il, jl, i, j

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ALLOCATE(JACOBI_DIAG(ib:il,jb:jl))

  DO j=jb,jl
    DO i=ib,il
      JACOBI_DIAG(i,j) = 1.0_precs/STENCIL%central(i,j)
    ENDDO
  ENDDO

END SUBROUTINE PREP_JACOBI

! Jacobi preconditioner
SUBROUTINE JACOBI(r)
  REAL(PREC), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, jb, il, jl, i, j

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  DO j=jb,jl
    DO i=ib,il
      r(i,j) = r(i,j)*JACOBI_DIAG(i,j)
    ENDDO
  ENDDO

END SUBROUTINE JACOBI

! Preparation of the ILU(0) preconditioner
SUBROUTINE PREP_ILU0(prototype)
  REAL(PREC), INTENT(in) :: prototype
  INTEGER :: ib, jb, il, jl, i, j
  REAL(PREC), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ALLOCATE(ILU0_DIAG(ib:il,jb:jl))

  ! Define some aliases
  uf => STENCIL%zonal
  vf => STENCIL%meridional
  ff => STENCIL%central

  ! Calculate diagonal elements of ILU(0) factorization
  ILU0_DIAG = 0.0_precs
  ILU0_DIAG(ib,jb) = ff(ib,jb)
  DO i=ib+1,il
    ILU0_DIAG(i,jb) = ff(i,jb) - uf(i-1,jb)**2/ILU0_DIAG(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    ILU0_DIAG(ib,j) = ff(ib,j) - vf(ib,j-1)**2/ILU0_DIAG(ib,j-1)
  ENDDO
  DO j=jb+1,jl
    DO i=ib+1,il
      ILU0_DIAG(i,j) = ff(i,j) - uf(i-1,j)**2/ILU0_DIAG(i-1,j) &
           - vf(i,j-1)**2/ILU0_DIAG(i,j-1)
    ENDDO
  ENDDO

END SUBROUTINE PREP_ILU0

! ILU(0) preconditioner
SUBROUTINE ILU0(r)
  REAL(PREC), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, jb, il, jl, i, j
  REAL(PREC), DIMENSION(:,:), POINTER :: uf, vf

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Define some aliases
  uf => STENCIL%zonal
  vf => STENCIL%meridional

  ! The overall system to solve is As = r <-> (LU)s = r

  ! 1.) Solve system Ly = r
  DO i=ib+1,il
    r(i,jb) = r(i,jb) + uf(i-1,jb)/ILU0_DIAG(i-1,jb)*r(i-1,jb)
  ENDDO
  DO j=jb+1,jl
    r(ib,j) = r(ib,j) + vf(ib,j-1)/ILU0_DIAG(ib,j-1)*r(ib,j-1)
  ENDDO
  DO j = jb+1,jl
    DO i = ib+1,il
      r(i,j) = r(i,j) + uf(i-1,j)/ILU0_DIAG(i-1,j)*r(i-1,j) &
           + vf(i,j-1)/ILU0_DIAG(i,j-1)*r(i,j-1)
    ENDDO
  ENDDO

  ! 2.) Solve system Us = y
  r(il,jl) = 1.0_precs/ILU0_DIAG(il,jl)*r(il,jl)
  DO i=il-1,ib,-1
    r(i,jl) = 1.0_precs/ILU0_DIAG(i,jl)*( r(i,jl) + uf(i,jl)*r(i+1,jl) )
  ENDDO
  DO j=jl-1,jb,-1
    r(il,j) = 1.0_precs/ILU0_DIAG(il,j)*( r(il,j) + vf(il,j)*r(il,j+1) )
  ENDDO
  DO j = jl-1,jb,-1
    DO i = il-1,ib,-1
      r(i,j) = 1.0_precs/ILU0_DIAG(i,j)*( r(i,j) + uf(i,j)*r(i+1,j)   &
           + vf(i,j)*r(i,j+1) )
    ENDDO
  ENDDO

END SUBROUTINE ILU0

! Symmetric SOR preconditioner
SUBROUTINE SSOR(r)

  REAL(PREC), INTENT(INOUT) :: r(:,:)
  INTEGER :: i, j, ib, jb, il, jl
  REAL(PREC):: ssorpar
  REAL(PREC), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! Get SSOR parameter from configuration
  ssorpar = REAL(get_ssor_param(), PREC)

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Define some aliases
  uf => STENCIL%zonal
  vf => STENCIL%meridional
  ff => STENCIL%central

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

END SUBROUTINE SSOR

! Calculates Matrix for Incomplete Cholesky Decomposition (ICC) with fillin p
! This routine makes use of a transformation from the barotropic stencil to a
! real matrix A with indices m_i, m_j. Variables s_i, s_j refer to the stencil
! coordinates, i, j are used according the actual context.
! The result L of L*L' ~= A is saved in ICC_C, ICC_W, ICC_S
SUBROUTINE PREP_ICC(p,  prototype, modify_opt)
  INTEGER, INTENT(IN) :: p  ! FILL-IN
  REAL(PREC), INTENT(IN) :: prototype
  LOGICAL, INTENT(IN), OPTIONAL :: modify_opt
  LOGICAL :: modify
  INTEGER :: ib, il, jb, jl, m_i, m_j, m_n, m_m, k
  REAL(PREC) :: tmp, d, cs
  REAL(PREC), DIMENSION(:,:), POINTER :: uf, vf, ff

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = KIND(prototype)

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Define some aliases
  uf => STENCIL%zonal
  vf => STENCIL%meridional
  ff => STENCIL%central

  ! Check and set optional arguments
  ! Check if modified ICC should be calculated
  IF ( PRESENT(modify_opt) ) THEN
    modify = modify_opt
  ELSE
    modify = .FALSE.
  ENDIF

  ! MAX(1,p-1) handles structurel exception for p=1
  ! ICC_C = Center, _W = West, _S = South in stencil
  ALLOCATE(ICC_C(ib:il,jb:jl), ICC_W(MAX(1,p-1),ib:il,jb:jl), &
       ICC_S(1+p,ib:il,jb:jl))

  ! Initialize
  ICC_C = 1.0_precs
  ICC_W = 0.0_precs
  ICC_S = 0.0_precs

  m_n = extent_size(STENCIL%extent)    ! dim of barotropic Matrix, m_ = Matrix
  m_m = extent_size(STENCIL%extent(1)) ! m_m is dist between diag and outer diag

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

  ! Calculate inverse of ICC_C to avoid division afterwards
  ICC_C = 1.0_precs/ICC_C

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
    REAL(PREC) :: res

    res = 0.0_precs

    ! return 0 if indices not in matrix (used if modify = True)
    IF ( i < 1 .OR. j < 1 .OR. i > m_n .OR. j > m_n ) RETURN

    ! diagonal
    IF ( i == j) THEN
      CALL coord_by_row(i, s_i, s_j)
      res = ICC_C(s_i, s_j)
      RETURN
    ENDIF

    ! secondary diagonal
    d = i - j
    IF ( d >= 1 .AND. d <= MAX(1,p-1) ) THEN
      CALL coord_by_row(i, s_i, s_j)
      res = ICC_W(d, s_i, s_j)
      RETURN
    ENDIF

    ! outer diagonal
    IF ( i - m_m > 0 ) THEN
      d = j - (i - m_m) + 1
      IF ( d >= 1 .AND. d <= p+1 ) THEN
        CALL coord_by_row(i, s_i, s_j)
        res = ICC_S(d, s_i, s_j)
        RETURN
      ENDIF
    ENDIF

  END FUNCTION get_value_of_L

  ! Set value in matrix L
  SUBROUTINE set_value_of_L(i,j,v)
    INTEGER, INTENT(IN) :: i, j
    REAL(PREC), INTENT(IN) :: v
    INTEGER :: d, s_i, s_j

    ! diagonal
    IF ( i == j) THEN
      CALL coord_by_row(i, s_i, s_j)
      ICC_C(s_i, s_j) = v
      RETURN
    ENDIF

    ! secondary diagonal
    d = i - j
    IF ( d >= 1 .AND. d <= MAX(1,p-1) ) THEN
      CALL coord_by_row(i, s_i, s_j)
      ICC_W(d, s_i, s_j) = v
      RETURN
    ENDIF

    ! outer diagonal
    IF ( i - m_m > 0 ) THEN
      d = j - (i - m_m) + 1
      IF ( d >= 1 .AND. d <= p+1 ) THEN
        CALL coord_by_row(i, s_i, s_j)
        ICC_S(d, s_i, s_j) = v
        RETURN
      ENDIF
    ENDIF

  END SUBROUTINE set_value_of_L

  ! Get value of virtual barotropic matrix A
  FUNCTION get_value_of_A(i,j) RESULT(res)
    INTEGER, INTENT(IN) :: i, j
    INTEGER :: s_i, s_j
    REAL(PREC) :: res

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

END SUBROUTINE PREP_ICC

! ICC(p) preconditioner
! Attention: This function silently assumes that r is an ordered real
! (not NaN nor +-Inf) outside of r(ib:il,jb:jl)
SUBROUTINE ICC(r)
  REAL(PREC), INTENT(INOUT) :: r(:,:)
  INTEGER :: ib, il, jb, jl, i, j, k, p, s_l, w_l

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Set halos to zero would be more secure at this point but affects
  ! the performance negatively. This should be done outside of this.
  !CALL clear_halos(r, stencil%extent)

  ! get the fill-in p of ICC
  p = SIZE(ICC_S,1) - 1
  s_l = p + 1 ! length of south diagonal
  w_l = MAX(1,p-1) ! length of west diagonal

  ! The overall system to solve is As = r <-> (LL')s = r

  ! 1.) Solve system Ly = r

  ! first column
  r(ib,jb) = r(ib,jb)*ICC_C(ib,jb)
  DO i=ib+1,il
    r(i,jb) = ( r(i,jb) - ICC_W(1,i,jb)*r(i-1,jb) )*ICC_C(i,jb)
  ENDDO

  ! all other columns
  DO j = jb+1,jl
    DO i = ib,il
      DO k = 1,s_l
        r(i,j) = r(i,j) - ICC_S(k,i,j)*r(i+k-1,j-1)
      ENDDO
      DO k = 1,w_l
        r(i,j) = r(i,j) - ICC_W(k,i,j)*r(i-k,j)
      ENDDO
      r(i,j) = r(i,j)*ICC_C(i,j)
    ENDDO
  ENDDO

  ! 2.) Solve system L's = y
  ! Backward substitution is done on a column level instead of row-wise

  ! all but first columns
  DO j = jl,jb+1,-1
    DO i = il,ib,-1
      r(i,j) = r(i,j)*ICC_C(i,j)
      DO k = w_l,1,-1
        r(i-k,j) = r(i-k,j) - r(i,j)*ICC_W(k,i,j) ! saxpy
      ENDDO
      DO k = s_l,1,-1
        r(i+k-1,j-1) = r(i+k-1,j-1) - r(i,j)*ICC_S(k,i,j) ! saxpy
      ENDDO
    ENDDO
  ENDDO

  ! first column
  DO i = il,ib+1,-1
    r(i,jb) = r(i,jb)*ICC_C(i,jb)
    r(i-1,jb) = r(i-1,jb) - r(i,jb)*ICC_W(1,i,jb) ! saxpy
  ENDDO
  r(ib,jb) = r(ib,jb)*ICC_C(ib,jb)

END SUBROUTINE ICC

! Prepare MICC(p) preconditioner
! Modified Incomplete Cholesky with fill-in p. The modification is that the
! diagonal entries of ICC(p) get modified in the way that the row sums of
! L*L' equal the row sums of A.
! MICC(p) only differs from ICC(p) in the preparation step. The application
! of MICC(p) is done by calling the subroutine icc.
SUBROUTINE PREP_MICC(p, prototype)
  INTEGER, INTENT(IN) :: p  ! FILL-IN
  REAL(PREC), INTENT(IN) :: prototype
  CALL prep_icc(p, prototype, .TRUE.)
END SUBROUTINE PREP_MICC

! Check if a certain preconditioner is already prepared
FUNCTION PRECOND_PREPARED(preconditioner, prototype) RESULT(ans)
  INTEGER, INTENT(IN) :: preconditioner
  REAL(PREC), INTENT(IN) :: prototype
  LOGICAL :: ans
  CHARACTER(len=25), PARAMETER :: msg_prefix = 'Selected preconditioner ['
  CHARACTER(len=17), PARAMETER :: msg_suffix = '] does not exist!'
  CHARACTER(len=LEN(msg_prefix) + LEN(msg_suffix) + 11) :: msg

  ans = KIND(prototype) == PREC
  SELECT CASE (preconditioner)
  CASE (NONE_PRECOND)
    ans = .TRUE.
  CASE (JACOBI_PRECOND)
    ans = ALLOCATED(JACOBI_DIAG)
  CASE (ILU0_PRECOND)
    ans = ALLOCATED(ILU0_DIAG)
  CASE (SSOR_PRECOND)
    ans = .NOT. ieee_is_nan(config%ssor_param)
  CASE (ICC_PRECOND,MICC_PRECOND)
    ans = ALLOCATED(ICC_C) .AND. ALLOCATED(ICC_W) .AND. ALLOCATED(ICC_S)
  CASE DEFAULT
    WRITE (msg, '(a,i0,a)') msg_prefix, preconditioner, msg_suffix
    CALL abort_ppm(msg, filename, __LINE__)
  END SELECT

END FUNCTION PRECOND_PREPARED

! Jacobi preconditioned matrix stencil operation
SUBROUTINE JACOBI_PRECOND_STENCIL(field, res_field)

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)

  ! Apply stencil first
  IF (.NOT. precond_prepared(JACOBI_PRECOND, 0._precs)) THEN
    CALL prep_jacobi(0._precs)
  ENDIF
  CALL apply_stencil(field, res_field)

  ! Apply Jacobi preconditioner
  CALL jacobi(res_field)

END SUBROUTINE JACOBI_PRECOND_STENCIL

! Jacobi preconditioned matrix stencil operation shifted by a value
SUBROUTINE JACOBI_PRECOND_SHIFTED_STENCIL(field, res_field)

  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

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
      res_field(i,j) = STENCIL%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE JACOBI_PRECOND_SHIFTED_STENCIL

! ILU(0) preconditioned matrix stencil operation
SUBROUTINE ILU0_PRECOND_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ILU0 preconditioner
  IF (.NOT. precond_prepared(ILU0_PRECOND, 0._precs)) THEN
    CALL prep_ilu0(0._precs)
  ENDIF
  CALL ilu0(res_field)

END SUBROUTINE ILU0_PRECOND_STENCIL

! ILU(0) preconditioned matrix stencil operation shifted by a value
SUBROUTINE ILU0_PRECOND_SHIFTED_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ILU0 preconditioner
  IF (.NOT. precond_prepared(ILU0_PRECOND, 0._precs)) THEN
    CALL prep_ilu0(0._precs)
  ENDIF
  CALL ilu0(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = STENCIL%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE ILU0_PRECOND_SHIFTED_STENCIL

! SSOR preconditioned matrix stencil operation
SUBROUTINE SSOR_PRECOND_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply SSOR preconditioner
  IF (.NOT. precond_prepared(SSOR_PRECOND, 0._precs)) THEN
    CALL abort_ppm("No SSOR parameter provided!", filename, __LINE__)
  ENDIF
  CALL ssor(res_field)

END SUBROUTINE SSOR_PRECOND_STENCIL

! SSOR preconditioned matrix stencil operation shifted by a value
SUBROUTINE SSOR_PRECOND_SHIFTED_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply SSOR preconditioner
  IF (.NOT. precond_prepared(SSOR_PRECOND, 0._precs)) THEN
    CALL abort_ppm("No SSOR parameter provided!", filename, __LINE__)
  ENDIF
  CALL ssor(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = STENCIL%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE SSOR_PRECOND_SHIFTED_STENCIL

! ICC(p) preconditioned matrix stencil operation
SUBROUTINE ICC_PRECOND_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ICC preconditioner
  IF (.NOT. precond_prepared(ICC_PRECOND, 0._precs)) THEN
    CALL prep_icc(config%icc_param, 0._precs)
  ENDIF
  CALL icc(res_field)

END SUBROUTINE ICC_PRECOND_STENCIL

! ICC(p) preconditioned matrix stencil operation shifted by a value
SUBROUTINE ICC_PRECOND_SHIFTED_STENCIL(field, res_field)
  REAL(PREC), INTENT(IN) :: field(:,:)
  REAL(PREC), INTENT(OUT) :: res_field(:,:)
  INTEGER :: i, j , ib, jb, il, jl

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  ! Define some variables for the ranges of the fields
  ib = extent_start(STENCIL%extent(1))
  jb = extent_start(STENCIL%extent(2))
  il = extent_end(STENCIL%extent(1))
  jl = extent_end(STENCIL%extent(2))

  ! Apply stencil first
  CALL apply_stencil(field, res_field)

  ! Apply ICC preconditioner
  IF (.NOT. precond_prepared(ICC_PRECOND, 0._precs)) THEN
    CALL prep_icc(config%icc_param, 0._precs)
  ENDIF
  CALL icc(res_field)

  ! Perform the shift
  DO j = jb, jl
    DO i = ib, il
      res_field(i,j) = STENCIL%shift*field(i,j) - res_field(i,j)
    END DO
  END DO

END SUBROUTINE ICC_PRECOND_SHIFTED_STENCIL
#undef filename
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
