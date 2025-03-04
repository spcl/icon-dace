! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_fast_math_lib

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: vec_div
  PUBLIC :: vec_aint
  PUBLIC :: vec_anint
  PUBLIC :: vec_exp
  PUBLIC :: vec_expm1
  PUBLIC :: vec_log
  PUBLIC :: vec_logp1
  PUBLIC :: vec_log10
  PUBLIC :: vec_pow
  PUBLIC :: vec_sqrt
  PUBLIC :: vec_cbrt
  PUBLIC :: vec_qdrt
  PUBLIC :: vec_rec
  PUBLIC :: vec_rsqrt
  PUBLIC :: vec_rcbrt
  PUBLIC :: vec_rqdrt
  PUBLIC :: vec_sincos
  PUBLIC :: vec_cos
  PUBLIC :: vec_acos
  PUBLIC :: vec_cosh
  PUBLIC :: vec_sin
  PUBLIC :: vec_asin
  PUBLIC :: vec_sinh
  PUBLIC :: vec_tan
  PUBLIC :: vec_atan2
  PUBLIC :: vec_tanh

  REAL(dp), PARAMETER, PRIVATE :: onethird = 1.0_dp/3.0_dp
  REAL(dp), PARAMETER, PRIVATE :: onefourth = 0.25_dp

CONTAINS

  SUBROUTINE vec_div(x, y, z, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: y(:)
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: z(:)
    INTEGER :: vec_size
    vec_size = SIZE(z)
    IF (PRESENT(n)) vec_size = n


    z(1:vec_size) = x(1:vec_size)/y(1:vec_size)

  END SUBROUTINE vec_div

  SUBROUTINE vec_aint(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = aint(x(1:vec_size))

  END SUBROUTINE vec_aint

  SUBROUTINE vec_anint(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = anint(x(1:vec_size))

  END SUBROUTINE vec_anint

  SUBROUTINE vec_exp(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = exp(x(1:vec_size))

  END SUBROUTINE vec_exp

  SUBROUTINE vec_expm1(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = exp(x(1:vec_size))-1.0_dp

  END SUBROUTINE vec_expm1

  SUBROUTINE vec_log(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = log(x(1:vec_size))

  END SUBROUTINE vec_log

  SUBROUTINE vec_logp1(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = log(x(1:vec_size)+1)

  END SUBROUTINE vec_logp1

  SUBROUTINE vec_log10(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = log10(x(1:vec_size))

  END SUBROUTINE vec_log10

  SUBROUTINE vec_pow(x, y, z, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: y(:)
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: z(:)
    INTEGER :: vec_size
    vec_size = SIZE(z)
    IF (PRESENT(n)) vec_size = n

    z(1:vec_size) = x(1:vec_size)**y(1:vec_size)

  END SUBROUTINE vec_pow

  SUBROUTINE vec_sqrt(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = sqrt(x(1:vec_size))

  END SUBROUTINE vec_sqrt

  SUBROUTINE vec_cbrt(x, y, n, lopenacc)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    LOGICAL,  INTENT(in), OPTIONAL :: lopenacc
    INTEGER :: i, vec_size
    LOGICAL :: lzopenacc
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n


    y(1:vec_size) = x(1:vec_size)**onethird


  END SUBROUTINE vec_cbrt

  SUBROUTINE vec_qdrt(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = x(1:vec_size)**onefourth

  END SUBROUTINE vec_qdrt

  SUBROUTINE vec_rec(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = 1.0_dp/x(1:vec_size)

  END SUBROUTINE vec_rec

  SUBROUTINE vec_rsqrt(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = 1.0_dp/sqrt(x(1:vec_size))

  END SUBROUTINE vec_rsqrt

  SUBROUTINE vec_rcbrt(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = 1.0_dp/x(1:vec_size)**onethird

  END SUBROUTINE vec_rcbrt

  SUBROUTINE vec_rqdrt(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = 1.0_dp/x(1:vec_size)**onefourth

  END SUBROUTINE vec_rqdrt

  SUBROUTINE vec_sincos(x, ys, yc, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: ys(:)
    REAL(dp), INTENT(inout) :: yc(:)
    INTEGER :: vec_size
    vec_size = SIZE(ys)
    IF (PRESENT(n)) vec_size = n

    ys(1:vec_size) = sin(x(1:vec_size))
    yc(1:vec_size) = cos(x(1:vec_size))
    

  END SUBROUTINE vec_sincos

  SUBROUTINE vec_cos(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = cos(x(1:vec_size))

  END SUBROUTINE vec_cos

  SUBROUTINE vec_acos(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = acos(x(1:vec_size))

  END SUBROUTINE vec_acos

  SUBROUTINE vec_cosh(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = cosh(x(1:vec_size))

  END SUBROUTINE vec_cosh

  SUBROUTINE vec_sin(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = sin(x(1:vec_size))

  END SUBROUTINE vec_sin

  SUBROUTINE vec_asin(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = asin(x(1:vec_size))

  END SUBROUTINE vec_asin

  SUBROUTINE vec_sinh(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = sinh(x(1:vec_size))

  END SUBROUTINE vec_sinh

  SUBROUTINE vec_tan(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = tan(x(1:vec_size))

  END SUBROUTINE vec_tan

  SUBROUTINE vec_atan2(x, y, z, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: y(:)
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: z(:)
    INTEGER :: vec_size
    vec_size = SIZE(z)
    IF (PRESENT(n)) vec_size = n

    z(1:vec_size) = atan2(x(1:vec_size),y(1:vec_size))

  END SUBROUTINE vec_atan2

  SUBROUTINE vec_tanh(x, y, n)
    INTEGER,  INTENT(in), OPTIONAL :: n
    REAL(dp), INTENT(in)    :: x(:)
    REAL(dp), INTENT(inout) :: y(:)
    INTEGER :: vec_size
    vec_size = SIZE(y)
    IF (PRESENT(n)) vec_size = n

    y(1:vec_size) = tanh(x(1:vec_size))

  END SUBROUTINE vec_tanh

END MODULE mo_fast_math_lib
