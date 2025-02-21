!>
!! @file linear_algebra_multi.f90
!! @author Florian Wilhelm
!!
!! @copyright Copyright  (C)  2010  Florian Wilhelm <Florian.Wilhelm@kit.edu>
!!
!! @version 1.0
!! @author Florian Wilhelm <Florian.Wilhelm@kit.edu>
!
! Keywords: scales ppm solver linear algebra residual
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
FUNCTION ARR_DOTPRODCUT1(x, global_opt) RESULT(ans)
  REAL(PREC), INTENT(IN) :: x(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: global_opt
  REAL(PREC) :: ans
  LOGICAL :: global

  IF ( PRESENT(global_opt) ) THEN
    global = global_opt
  ELSE
    global = config%do_exchange
  ENDIF

  ans = SUM(x**2)

  IF (global) ans = global_sum(ans)

END FUNCTION ARR_DOTPRODCUT1

!> Array-wise dotproduct for two different matrices
FUNCTION ARR_DOTPRODCUT2(x, y, global_opt) RESULT(ans)



  REAL(PREC), INTENT(IN) :: x(:,:), y(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: global_opt
  REAL(PREC) :: ans
  REAL(PREC), ALLOCATABLE :: xy(:)
  LOGICAL :: global
  INTEGER :: j, je



  IF ( PRESENT(global_opt) ) THEN
    global = global_opt
  ELSE
    global = config%do_exchange
  ENDIF




  je = SIZE(x,2)


  ALLOCATE(xy(je))

  DO j=1,je
    xy(j) = DOT_PRODUCT(x(:,j), y(:,j))
  ENDDO

  ans = SUM(xy)
  DEALLOCATE(xy)



  IF (global) ans = global_sum(ans)

END FUNCTION ARR_DOTPRODCUT2

! Array-wise 2-norm
FUNCTION ARR_NORM_2(x, global_opt) RESULT(ans)
  REAL(PREC), INTENT(IN) :: x(:,:)
  LOGICAL, INTENT(IN), OPTIONAL :: global_opt
  REAL(PREC) :: ans
  LOGICAL :: global

  IF ( PRESENT(global_opt) ) THEN
    global = global_opt
  ELSE
    global = config%do_exchange
  ENDIF

  ans = SQRT(arr_dotproduct(x, global))

END FUNCTION ARR_NORM_2

! Sums a variable globally up
FUNCTION GLOBAL_SUM(summand, comm_opt) RESULT(all_sum)





  REAL(PREC), INTENT(in) :: summand
  INTEGER, OPTIONAL, INTENT(in) :: comm_opt
  REAL(PREC) :: all_sum
  all_sum = summand

END FUNCTION GLOBAL_SUM

! Calculate absolute residual of system
! If global_opt is .TRUE. then calculate the global absolute residual
FUNCTION CALC_ABS_RES(A, b, x, ext_x, global_opt) RESULT(abs_res)
  REAL(PREC), INTENT(IN) :: x(:,:)
  REAL(PREC), INTENT(IN) :: b(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)      ! extent of x
  LOGICAL, INTENT(IN), OPTIONAL :: global_opt
  REAL(PREC), ALLOCATABLE :: r(:,:), Ax(:,:)
  INTEGER :: ie, je, ib, jb, il, jl
  REAL(PREC) :: abs_res
  LOGICAL :: global

  INTERFACE
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
  END INTERFACE

  IF ( PRESENT(global_opt) ) THEN
    global = global_opt
  ELSE
    global = config%do_exchange
  ENDIF

  ie = SIZE(x,1)
  je = SIZE(x,2)
  ib = extent_start(ext_x(1))
  jb = extent_start(ext_x(2))
  il = extent_end(ext_x(1))
  jl = extent_end(ext_x(2))

  ALLOCATE(r(ie,je), Ax(ie,je))

  CALL A(x, Ax)
  r(ib:il,jb:jl) = b(ib:il,jb:jl) - Ax(ib:il,jb:jl)

  abs_res = arr_norm_2(r(ib:il,jb:jl), global)

  DEALLOCATE(Ax, r)

END FUNCTION CALC_ABS_RES

! Calculate relative residual of system
! If global_opt is .TRUE. then calculate the global relative residual
FUNCTION CALC_REL_RES(A, b, x, ext_x, global_opt, x0_opt) RESULT(rel_res)
  REAL(PREC), INTENT(IN) :: x(:,:)
  REAL(PREC), INTENT(IN) :: b(:,:)
  TYPE(extent), INTENT(IN) :: ext_x(:)      ! extent of x
  LOGICAL, INTENT(IN), OPTIONAL :: global_opt
  REAL(PREC), INTENT(IN), OPTIONAL :: x0_opt(:,:)
  REAL(PREC), ALLOCATABLE :: x0(:,:)
  INTEGER :: ie, je
  REAL(PREC) :: rel_res, numer, denom

  ! PREC suffix (precs) as constant to avoid Preprocessor problems
  INTEGER, PARAMETER :: precs = PREC

  INTERFACE
    SUBROUTINE A(field, res_field)
      USE ppm_std_type_kinds, ONLY: PREC
      REAL(PREC), INTENT(IN) :: field(:,:)
      REAL(PREC), INTENT(OUT) :: res_field(:,:)
    END SUBROUTINE A
  END INTERFACE

  IF ( PRESENT(x0_opt) ) THEN
    x0 = x0_opt
  ELSE
    ie = SIZE(x,1)
    je = SIZE(x,2)
    ALLOCATE(x0(ie,je))
    x0 = 0.0_precs
  ENDIF

  IF ( PRESENT(global_opt) ) THEN
    numer = calc_abs_res(A, b, x, ext_x, global_opt)
    denom = calc_abs_res(A, b, x0, ext_x, global_opt)
  ELSE
    numer = calc_abs_res(A, b, x, ext_x)
    denom = calc_abs_res(A, b, x0, ext_x)
  ENDIF

  rel_res = numer/denom

  IF (.NOT. PRESENT(x0_opt)) THEN
    DEALLOCATE(x0)
  ENDIF

END FUNCTION CALC_REL_RES
!
! Local Variables:
! mode: f90
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
