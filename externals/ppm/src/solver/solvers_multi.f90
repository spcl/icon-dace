!>
!! @file solvers_multi.f90
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
! Preconditioned CG-Method with A*x=b and startvalue x
FUNCTION PRECOND_CG_METHOD(A, b, x, ext_x, precond, exchange, tol_opt, &
     maxiter_opt) RESULT(kiter)
  REAL(PREC), INTENT(IN) :: b(:,:)                ! right-hand-side b
  REAL(PREC), INTENT(INOUT) :: x(:,:)             ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(PREC), OPTIONAL, INTENT(IN) :: tol_opt     ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
  INTEGER :: jb, ib, il, jl, ie, je, kiter, maxiter



  REAL(PREC) :: delta, delta0, delta_old, dAd, alpha, beta, tol

  ! Auxiliary fields for CG-Method
  REAL(PREC), ALLOCATABLE :: d(:,:), r(:,:), Ad(:,:), s(:,:)

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Preconditioner
    SUBROUTINE precond(r)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: r(:,:)
    END SUBROUTINE precond
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  ! Check and set optional arguments
  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt**2 ! **2, because we compare later the square of a norm
  ELSE
    tol = REAL(config%tol**2, PREC)
  ENDIF

  ! Define some variables for the ranges of the fields
  ie = SIZE(x,1)
  je = SIZE(x,2)
  ib = extent_start(ext_x(1))
  jb = extent_start(ext_x(2))
  il = extent_end(ext_x(1))
  jl = extent_end(ext_x(2))




  ALLOCATE(d(ie,je), r(ie,je), Ad(ie,je), s(ie,je))

  ! Initialize boundary/halos of d and s
  CALL clear_halos(d, ext_x)
  CALL clear_halos(s, ext_x)

  ! Initialise residual and direction
  CALL A(x, Ad)

  r(ib:il,jb:jl) = b(ib:il,jb:jl) - Ad(ib:il,jb:jl)




  ! Apply Preconditioner
  d(ib:il,jb:jl) = r(ib:il,jb:jl)
  CALL precond(d)

  ! Calculate scalar product of residual and direction
  delta = arr_dotproduct(r(ib:il,jb:jl), d(ib:il,jb:jl))

  delta0 = delta
  kiter = 0

  DO WHILE (kiter < maxiter .AND. delta > tol*delta0)
    IF (config%do_exchange) CALL exchange(d,'precond_cg_method 1')

    kiter = kiter + 1

    ! Initialise auxiliary variables
    CALL A(d, Ad)

    dAd = arr_dotproduct( d(ib:il,jb:jl), Ad(ib:il,jb:jl) )

    alpha = delta/dAd

    ! One step of preconditioned CG-method

    x(ib:il,jb:jl) = x(ib:il,jb:jl) + alpha*d(ib:il,jb:jl)
    r(ib:il,jb:jl) = r(ib:il,jb:jl) - alpha*Ad(ib:il,jb:jl)





    ! Apply Preconditioner
    s(ib:il,jb:jl) = r(ib:il,jb:jl)
    CALL precond(s)

    ! Calculate scalar product of residual and direction
    delta_old = delta
    delta = arr_dotproduct( r(ib:il,jb:jl), s(ib:il,jb:jl) )

    beta = delta/delta_old

    !write(*,*) 'PCG iteration: ', kiter, ' current delta: ', delta

    ! Compute new direction d

    d(ib:il,jb:jl) = s(ib:il,jb:jl) + beta*d(ib:il,jb:jl)






  ENDDO

  DEALLOCATE(d, r, s, Ad)

  ! Communicate x (only d was used the whole time)
  IF (config%do_exchange) CALL exchange(x, 'precond_cg_method 2')

END FUNCTION PRECOND_CG_METHOD

! CG-Method with A*x=b and startvalue x
FUNCTION CG_METHOD(A, b, x, ext_x, exchange, tol_opt, maxiter_opt) RESULT(kiter)

  REAL(PREC), INTENT(IN) :: b(:,:)                  ! right-hand-side b
  REAL(PREC), INTENT(INOUT) :: x(:,:)               ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(PREC), OPTIONAL, INTENT(IN) :: tol_opt       ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
  INTEGER :: kiter                                ! number of iterations

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  IF ( PRESENT(tol_opt) ) THEN
    IF ( PRESENT(maxiter_opt) ) THEN
      kiter = precond_cg_method(A, b, x, ext_x, IDENTITY, exchange, tol_opt, &
           maxiter_opt)
    ELSE
      kiter = precond_cg_method(A, b, x, ext_x, IDENTITY, exchange, tol_opt)
    ENDIF
  ELSE
    kiter = precond_cg_method(A, b, x, ext_x, IDENTITY, exchange)
  ENDIF

END FUNCTION CG_METHOD

! Additive Schwarz method which uses a local solver for each partition
FUNCTION SCHWARZ_METHOD(A, b, x, ext_x, local_solver, exchange, tol_opt, &
     maxiter_opt, omega_opt) RESULT(kiter)
  REAL(PREC), INTENT(IN) :: b(:,:)                  ! right-hand-side b
  REAL(PREC), INTENT(INOUT) :: x(:,:)               ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)              ! extent of x
  REAL(PREC), OPTIONAL, INTENT(IN) :: tol_opt       ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt      ! maximum iterations
  REAL(PREC), OPTIONAL, INTENT(IN) :: omega_opt     ! relaxation parameter
  INTEGER :: kiter, maxiter, iiter, ie, je, jb, ib, jl, il!, ierror



  !INTEGER :: my_rank
  REAL(PREC) :: tol, numer, denom, omega!, local_relres
  REAL(PREC), ALLOCATABLE :: r(:,:), dx(:,:), Ax(:,:)

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    !> the solver used by Schwarz method
    FUNCTION local_solver(A, b, x, ext_x, exchange, tol_opt, maxiter_opt) &
         RESULT(kiter)
      USE ppm_std_type_kinds, ONLY: PREC
      USE ppm_extents, ONLY: extent

      REAL(PREC), INTENT(IN) :: b(:,:)                ! right-hand-side b
      REAL(PREC), INTENT(INOUT) :: x(:,:)             ! startvalue and result
      TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
      REAL(PREC), OPTIONAL, INTENT(IN) :: tol_opt     ! tolerance for residual
      INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
      INTEGER :: kiter                                ! needed iterations
      INTERFACE
        ! Matrix-Vector multiplication of the linear system to solve
        SUBROUTINE A(field, res_field)
          USE ppm_std_type_kinds, ONLY: PREC
          REAL(PREC), INTENT(IN) :: field(:,:)
          REAL(PREC), INTENT(OUT) :: res_field(:,:)
        END SUBROUTINE A
        ! Function to exchange boundaries if necessary
        SUBROUTINE exchange(a0, text)
          USE ppm_std_type_kinds, ONLY: PREC
          REAL(PREC), INTENT(INOUT) :: a0(:,:)
          CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
        END SUBROUTINE exchange
      END INTERFACE
    END FUNCTION local_solver
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%schwarz_maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt
  ELSE
    tol = REAL(config%schwarz_tol, PREC)
  ENDIF
  IF ( PRESENT(omega_opt) ) THEN
    omega = omega_opt
  ELSE
    omega = REAL(config%schwarz_relax, PREC)
  ENDIF

  kiter = 0
  ie = SIZE(b,1)
  je = SIZE(b,2)
  ib = extent_start(ext_x(1))
  jb = extent_start(ext_x(2))
  il = extent_end(ext_x(1))
  jl = extent_end(ext_x(2))




  ALLOCATE(r(ie,je), dx(ie,je), Ax(ie,je))

  ! Initialise residual
  CALL A(x, Ax)

  r(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)




  denom = arr_norm_2(r(ib:il,jb:jl), .TRUE.)
  numer = denom

  DO WHILE (kiter < maxiter .AND. numer > tol*denom)
    kiter = kiter + 1

    dx = 0.0_precs

    ! Solve local problem according to solver as configured
    iiter = local_solver(A, r, dx, ext_x, exchange, REAL(config%tol, PREC), &
         config%maxiter)


    x(ib:il,jb:jl) = x(ib:il,jb:jl) + omega*dx(ib:il,jb:jl)




    CALL exchange(x, 'schwarz_method 1')

    CALL A(x, Ax)

    r(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)




!~     local_relres = calc_rel_res(A, b, x, ext_x, .FALSE.)
!~     WRITE(*,*) "CPU: ", my_rank, "kiter: ", kiter, "iiter: ", iiter, &
!~         "local_relres: ", local_relres
!~     IF( my_rank == 0 ) WRITE(*,*) "Global relres", kiter, ":", numer/denom

    IF( kiter >= config%schwarz_miniter .AND. &
         MOD(kiter, config%schwarz_checkrate) == 0 ) THEN
      numer = arr_norm_2(r(ib:il,jb:jl), .TRUE.)
    ENDIF
  END DO

  DEALLOCATE(r, dx, Ax)

END FUNCTION SCHWARZ_METHOD

! Non-preconditioned Chebyshev method
FUNCTION CHEBYSHEV_METHOD(A, b, x, ext_x, lambda_min, lambda_max, exchange, &
     tol_opt, maxiter_opt) RESULT(kiter)
  REAL(PREC), INTENT(IN) :: b(:,:)                ! right-hand-side b
  REAL(PREC), INTENT(INOUT) :: x(:,:)             ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)            ! extent of x
  REAL(PREC), INTENT(IN) :: lambda_min, lambda_max! smallest, largest abs. EV
  REAL(PREC), OPTIONAL, INTENT(IN) :: tol_opt     ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt    ! maximum iterations
  INTEGER :: kiter                                ! number of iterations

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  IF ( PRESENT(tol_opt) ) THEN
    IF ( PRESENT(maxiter_opt) ) THEN
      kiter = precond_chebyshev_method(A, b, x, ext_x, lambda_min, lambda_max, &
           IDENTITY, exchange, tol_opt, maxiter_opt)
    ELSE
      kiter = precond_chebyshev_method(A, b, x, ext_x, lambda_min, lambda_max, &
           IDENTITY, exchange, tol_opt)
    ENDIF
  ELSE
    kiter = precond_chebyshev_method(A, b, x, ext_x, lambda_min, lambda_max, &
         IDENTITY, exchange)
  ENDIF

END FUNCTION CHEBYSHEV_METHOD

! Preconditioned Chebyshev acceleration according to Saad "Iterative Methods
! for Sparse Linear Systems", Algorithm 12.1
FUNCTION PRECOND_CHEBYSHEV_METHOD(A, b, x, ext_x, lambda_min, lambda_max, &
     precond, exchange, tol_opt, maxiter_opt) RESULT(kiter)
  REAL(PREC), INTENT(IN) :: b(:,:)                  ! right-hand-side b
  REAL(PREC), INTENT(INOUT) :: x(:,:)               ! startvalue and result
  TYPE(extent), INTENT(IN) :: ext_x(:)              ! extent of x
  REAL(PREC), INTENT(IN) :: lambda_min, lambda_max  ! smallest, largest abs. EV
  REAL(PREC), OPTIONAL, INTENT(IN) :: tol_opt       ! tolerance for residual
  INTEGER, OPTIONAL, INTENT(IN) :: maxiter_opt      ! maximum iterations
  INTEGER :: jb, ib, il, jl, ie, je, kiter, maxiter



  REAL(PREC) :: rho_old, rho_new, theta, delta, sigma, tol, denom, numer
  ! Auxiliary fields for CG-Method
  REAL(PREC), ALLOCATABLE :: r(:,:), Ax(:,:), d(:,:)

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  INTERFACE
    ! Matrix-Vector multiplication of the linear system to solve
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
    ! Preconditioner
    SUBROUTINE precond(r)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: r(:,:)
    END SUBROUTINE precond
    ! Function to exchange boundaries if necessary
    SUBROUTINE exchange(a0, text)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(INOUT) :: a0(:,:)
      CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: text
    END SUBROUTINE exchange
  END INTERFACE

  ! Check and set optional arguments
  IF ( PRESENT(maxiter_opt) ) THEN
    maxiter = maxiter_opt
  ELSE
    maxiter = config%maxiter
  ENDIF
  IF ( PRESENT(tol_opt) ) THEN
    tol = tol_opt**2 ! **2, because we compare later the square of a norm
  ELSE
    tol = REAL(config%tol**2, PREC)
  ENDIF

  ie = SIZE(x,1)
  je = SIZE(x,2)
  kiter = 0
  ib = extent_start(ext_x(1))
  jb = extent_start(ext_x(2))
  il = extent_end(ext_x(1))
  jl = extent_end(ext_x(2))




  theta = (lambda_max + lambda_min) / 2.0_precs
  delta = (lambda_max - lambda_min) / 2.0_precs

  ! Check for trivial solution
  ! use workaround because gfortran warns about comparing real with ==

  IF (delta >= 0.0_precs .AND. delta <= 0.0_precs) THEN
    IF (lambda_max  >= 0.0_precs .AND. lambda_max <= 0.0_precs) THEN




      x(ib:il,jb:jl) = 0.0_precs
    ELSE
      x(ib:il,jb:jl) = b(ib:il,jb:jl)/lambda_max
    ENDIF
    RETURN
  ENDIF

  ALLOCATE(r(ie,je), Ax(ie,je), d(ie,je))

  ! Initialize boundary/halos of r
  CALL clear_halos(r, ext_x)

  ! Residual and direction for first iteration
  CALL A(x, Ax)

  r(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)



  CALL precond(r)

  d(ib:il,jb:jl) = 1.0_precs/theta*r(ib:il,jb:jl)




  ! Denominator for relativ residual
  denom = arr_dotproduct(r(ib:il,jb:jl))
  numer = denom

  sigma = theta/delta
  rho_old = 1.0_precs/sigma

  DO WHILE (kiter < maxiter .AND. numer > tol*denom)
    kiter = kiter + 1

    ! x_new = x_old + d_old

    x(ib:il,jb:jl) = x(ib:il,jb:jl) + d(ib:il,jb:jl)



    IF (config%do_exchange) CALL exchange(x, 'chebyshev 1')

    ! r_new = r_old - A*d_old = b - A*x_new
    CALL A(x, Ax)

    r(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl) ! = r - Ad in Saad



    CALL precond(r)

    ! update rho
    rho_new = 1.0_precs/(2.0_precs*sigma - rho_old)

    ! d_new = rho_new*rho_old * d_old + (2*rho_new)/delta * r_new

    d(ib:il,jb:jl) = rho_new*rho_old*d(ib:il,jb:jl) &
         + (2.0_precs * rho_old)/delta*r(ib:il,jb:jl)





    rho_old = rho_new

    ! calculate current residual
    IF( kiter >= config%cheby_miniter .AND. &
         MOD(kiter, config%cheby_checkrate) == 0 ) THEN
      numer = arr_dotproduct(r(ib:il,jb:jl))
    ENDIF

    !WRITE(*,*) 'CHEBY: kiter, ', kiter, 'rel_res', numer/denom
  ENDDO

  DEALLOCATE(r, Ax, d)

END FUNCTION PRECOND_CHEBYSHEV_METHOD
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
