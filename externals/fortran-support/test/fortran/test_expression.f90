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
    REAL(wp), TARGET :: val
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
      CALL ASSERT_EQUAL(.TRUE., ALL(SIN(z(:, 1, :)*z_sfc) == val_3D(:, 1, :)))
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
