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

MODULE test_mo_expression
  USE FORTUTF
  USE mo_expression
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64

CONTAINS

  SUBROUTINE TEST_expression_simple
    TYPE(expression) :: formula
    REAL(wp), TARGET :: val, x
    REAL(wp), POINTER :: ptr_val
    REAL(wp) :: ref

    ptr_val => val

    CALL TAG_TEST("TEST_div_mult")
    formula = expression("20*10/20")
    ref = 20._wp*10._wp/20._wp
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_ALMOST_EQUAL(val, ref)

    CALL TAG_TEST("TEST_sin_cos")
    formula = expression("sin(10.)*cos(5.)*1.2")
    CALL formula%evaluate(ptr_val)
    ref = SIN(10._wp)*COS(5._wp)*1.2_wp
    CALL ASSERT_ALMOST_EQUAL(val, ref)

    CALL TAG_TEST("TEST_log_pow")
    formula = expression("log(8)*10^5")
    ref = LOG(8._wp)*10._wp**5
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_ALMOST_EQUAL(val, ref)

    CALL TAG_TEST("TEST_plus_minus")
    formula = expression("5 - 1 + 3")
    ref = REAL(5 - 1 + 3, wp)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_ALMOST_EQUAL(val, ref)

    CALL TAG_TEST("TEST_exp_erf")
    formula = expression("exp(1.) + erf(0.76)")
    ref = EXP(1._wp) + ERF(0.76_wp)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_ALMOST_EQUAL(val, ref)

    CALL TAG_TEST("TEST_max")
    formula = expression("max(1., 0.)")
    ref = 1._wp
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(val, ref)

    CALL TAG_TEST("TEST_min")
    formula = expression("min(1., 0.)")
    ref = 0._wp
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(val, ref)

    CALL TAG_TEST("TEST_gt")
    formula = expression("if(1 > 0, 1, 0)")
    ref = 1._wp
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(val, ref)

    CALL TAG_TEST("TEST_lt")
    formula = expression("if(1 < 0, 1, 0)")
    ref = 0._wp
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(val, ref)

    CALL TAG_TEST("TEST_sqrt")
    formula = expression("sqrt(2.) + sqrt(3.)")
    ref = SQRT(2._wp) + SQRT(3._wp)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_ALMOST_EQUAL(val, ref)

    CALL TAG_TEST("TEST_link")
    CALL RANDOM_NUMBER(x)
    formula = expression("[x]")
    CALL formula%link("x", x)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_ALMOST_EQUAL(val, x)

  END SUBROUTINE

  SUBROUTINE TEST_expression_2D()

    REAL(wp), PARAMETER :: TOL = 5.e-6
    TYPE(expression) :: formula
    REAL(wp), TARGET :: val(5, 5)
    REAL(wp), POINTER :: ptr_val(:, :)
    REAL(wp), TARGET :: x(5, 5), y(5, 5), z(5, 5)
    REAL(wp) :: ref(5, 5)
    LOGICAL :: logl(5, 5)

    ptr_val => val

    CALL RANDOM_NUMBER(x)
    CALL RANDOM_NUMBER(y)
    CALL RANDOM_NUMBER(z)

    CALL TAG_TEST("TEST_plus_minus_2D")
    formula = expression("[x] - [y] + [z]")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    CALL formula%link("z", z)
    ref = x - y + z
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_arithmetic_2D")
    formula = expression("[x] * 10 + [y] / 5 - [z]")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    CALL formula%link("z", z)
    ref = x*10 + y/5 - z
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_log_pow_2D")
    formula = expression("log([x]) ^ 2")
    CALL formula%link("x", x)
    ref = LOG(x)**2
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_exp_erf_2D")
    formula = expression("exp([x]) + erf([y])")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = EXP(x) + ERF(y)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_sin_cos_2D")
    formula = expression("sin([x]) + cos([y])")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = SIN(x) + COS(y)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_sqrt_2D")
    formula = expression("sqrt([x])")
    CALL formula%link("x", x)
    ref = SQRT(x)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_max_2D")
    formula = expression("max([x] - 0.5, 0)")
    CALL formula%link("x", x)
    ref = MAX(x - 0.5_wp, 0.0_wp)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_min_2D")
    formula = expression("min([x], [y])")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = MIN(x, y)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_gt_2D")
    formula = expression("[x] > 0.5")
    CALL formula%link("x", x)
    ref = REAL(ABS(FLOOR(0.5_wp - x)), KIND(REAL(wp)))
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(0.0_wp, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_lt_2D")
    formula = expression("[x] < 0.3")
    CALL formula%link("x", x)
    ref = REAL(1 + FLOOR(0.3_wp - x), KIND(REAL(wp)))
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(0.0_wp, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_if_2D")
    formula = expression("if([x] > [y], 1.0, 0.0)")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = REAL(ABS(FLOOR(y - x)), KIND(REAL(wp)))
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(0.0_wp, MAXVAL(ABS(val - ref)))

    CALL formula%finalize()

  END SUBROUTINE

  SUBROUTINE TEST_expression_3D()

    REAL(wp), PARAMETER :: TOL = 5.e-6
    TYPE(expression) :: formula
    REAL(wp), TARGET :: val(3, 3, 3)
    REAL(wp), POINTER :: ptr_val(:, :, :)
    REAL(wp), TARGET :: x(3, 3, 3), y(3, 3, 3), z(3, 3, 3)
    REAL(wp) :: ref(3, 3, 3)

    ptr_val => val

    CALL RANDOM_NUMBER(x)
    CALL RANDOM_NUMBER(y)
    CALL RANDOM_NUMBER(z)

    CALL TAG_TEST("TEST_plus_minus_3D")
    formula = expression("[x] - [y] + [z]")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    CALL formula%link("z", z)
    ref = x - y + z
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_arithmetic_3D")
    formula = expression("[x] * 10 + [y] / 5 - [z]")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    CALL formula%link("z", z)
    ref = x*10 + y/5 - z
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_log_pow_3D")
    formula = expression("log([x]) ^ 2")
    CALL formula%link("x", x)
    ref = LOG(x)**2
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_exp_erf_3D")
    formula = expression("exp([x]) + erf([y])")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = EXP(x) + ERF(y)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_sin_cos_3D")
    formula = expression("sin([x]) + cos([y])")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = SIN(x) + COS(y)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_sqrt_3D")
    formula = expression("sqrt([x])")
    CALL formula%link("x", x)
    ref = SQRT(x)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_max_3D")
    formula = expression("max([x] - 0.5, 0)")
    CALL formula%link("x", x)
    ref = MAX(x - 0.5_wp, 0.0_wp)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_min_3D")
    formula = expression("min([x], [y])")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = MIN(x, y)
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_gt_3D")
    formula = expression("[x] > 0.5")
    CALL formula%link("x", x)
    ref = REAL(ABS(FLOOR(0.5_wp - x)), KIND(REAL(wp)))
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(0.0_wp, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_lt_3D")
    formula = expression("[x] < 0.3")
    CALL formula%link("x", x)
    ref = REAL(1 + FLOOR(0.3_wp - x), KIND(REAL(wp)))
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(0.0_wp, MAXVAL(ABS(val - ref)))

    CALL TAG_TEST("TEST_if_3D")
    formula = expression("if([x] > [y], 1.0, 0.0)")
    CALL formula%link("x", x)
    CALL formula%link("y", y)
    ref = REAL(ABS(FLOOR(y - x)), KIND(REAL(wp)))
    CALL formula%evaluate(ptr_val)
    CALL ASSERT_EQUAL(0.0_wp, MAXVAL(ABS(val - ref)))

    CALL formula%finalize()

  END SUBROUTINE

  SUBROUTINE TEST_expression_complex()

    REAL(wp), PARAMETER :: pi = 4._wp*ATAN(1._wp)
    REAL(wp), PARAMETER :: TOL = 5.e-6
    LOGICAL                  :: lassert
    INTEGER                  :: i
    TYPE(expression)         :: formula
    REAL(wp), POINTER :: val => NULL()
    REAL(wp), POINTER :: val_2D(:, :) => NULL()
    REAL(wp), POINTER :: val_3D(:, :, :) => NULL()
    REAL(wp), TARGET  :: z(2, 3, 2), z_sfc(2, 2)
    CHARACTER(len=100) :: label

    ! create some dummy data
    z_sfc(1, :) = [1, 2]
    z_sfc(2, :) = [3, 4]
    z(:, 1, :) = 1
    z(:, 2, :) = 2
    z(:, 3, :) = 3

    ! --- 3D array example:
    ! create a formula expression
    formula = expression("sin([z] * [z_sfc])")
    CALL formula%link("z", z)
    CALL formula%link("z_sfc", z_sfc)

    ! evaluate the expression:
    CALL formula%evaluate(val_3D)

    ! check results and clean-up
    DO i = 1, 3
      WRITE (label, '(A,I1)') 'TEST_exp_sin_z_sfc_', i
      CALL TAG_TEST(label)
      CALL ASSERT_EQUAL(.TRUE., ALL(SIN(z(:, i, :)*z_sfc) == val_3D(:, i, :)))
    END DO

    ! clean-up
    DEALLOCATE (val_3D)
    CALL formula%finalize()

    ! --- 2D example:
    formula = expression("if([z_sfc] > 2., [z_sfc], 0. )")
    CALL formula%link("z_sfc", z_sfc)
    CALL formula%evaluate(val_2D)
    CALL TAG_TEST('TEST_exr_z_sfc_z_sfc')
    CALL ASSERT_EQUAL(.TRUE., ALL(val_2d == MERGE(z_sfc, 0._wp, z_sfc > 2.)))

    ! clean-up
    DEALLOCATE (val_2D)
    CALL formula%finalize()

    ! --- scalar examples
    formula = expression("sin(45*pi/180.) * 10 + 5")
    CALL formula%evaluate(val)
    CALL TAG_TEST('TEST_exr_sin(45*pi/180.)')
    CALL ASSERT_GREATER_THAN(TOL, ABS(val - (SIN(45.*pi/180.)*10.+5.)))
    DEALLOCATE (val)
    CALL formula%finalize()

    formula = expression("3./2.*pi")
    CALL formula%evaluate(val)
    CALL TAG_TEST('TEST_exr_3./2.*pi')
    CALL ASSERT_GREATER_THAN(TOL, ABS(val - (3./2.*pi)))
    DEALLOCATE (val)
    CALL formula%finalize()

    formula = expression("sin(max(2,3)/3 * 3.1415)")
    CALL formula%evaluate(val)
    CALL TAG_TEST('TEST_exr_sin(max(2,3)/3 * 3.1415)')
    CALL ASSERT_GREATER_THAN(TOL, ABS(val - (SIN(MAX(2, 3)/3*3.1415))))
    DEALLOCATE (val)
    CALL formula%finalize()

    formula = expression("3 + 4*2/(1 - 5)^2^3")
    CALL formula%evaluate(val)
    CALL TAG_TEST('TEST_exr_3 + 4*2/(1 - 5)^2^3')
    CALL ASSERT_GREATER_THAN(TOL, ABS(val - (3.+8./4096.)))
    DEALLOCATE (val)
    CALL formula%finalize()

    ! --- scalar example, written to pre-allocated result:
    ALLOCATE (REAL(wp)::val_2D(2, 2))
    formula = expression("sin(45*pi/180.) * 10 + 5")
    CALL formula%evaluate(val_2D)
    CALL TAG_TEST('TEST_exr_sin(45*pi/180.) * 10 + 5')
    CALL ASSERT_GREATER_THAN(TOL, MAXVAL(ABS(val_2D - (SIN(45.*pi/180.)*10.+5.))))
    DEALLOCATE (val_2D)
    CALL formula%finalize()
  END SUBROUTINE TEST_expression_complex

END MODULE test_mo_expression
