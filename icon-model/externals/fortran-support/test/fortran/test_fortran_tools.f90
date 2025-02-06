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

MODULE test_fortran_tools
  USE FORTUTF
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: dp => real64, &
    &                                      sp => real32, &
    &                                      i4 => int32
  USE mo_fortran_tools
  USE mo_exception, ONLY: init_logger, finish
  USE helpers, ONLY: open_logfile, open_new_logfile, custom_exit_dummy

CONTAINS

  SUBROUTINE Test_assign_if_present_character
    CHARACTER(len=1) :: x, y

    x = 'x'
    y = 'y'
    CALL TAG_TEST("Test_assign_with_present_character")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(x, y)

    x = ' '
    CALL TAG_TEST("Test_assign_with_empty_character")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(x, ' ')

    CALL TAG_TEST("Test_assign_no_present_character")
    CALL assign_if_present(y)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_logical
    LOGICAL :: x, y

    x = .TRUE.
    y = .FALSE.
    CALL TAG_TEST("Test_assign_with_present_logical")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(x, y)

    CALL TAG_TEST("Test_assign_no_present_logical")
    CALL assign_if_present(y)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_logicals
    LOGICAL :: x(5) = .TRUE.
    LOGICAL :: y(5) = (/.FALSE., .TRUE., .FALSE., .TRUE., .FALSE./)

    CALL TAG_TEST("Test_assign_with_present_logicals")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(assert_logical_array(x, y), .TRUE.)

    CALL TAG_TEST("Test_assign_no_present_logicals")
    CALL assign_if_present(y)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_integer
    INTEGER :: x, y

    x = 1
    y = 2
    CALL TAG_TEST("Test_assign_with_present_integer")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(x, y)

    CALL TAG_TEST("Test_assign_no_present_integer")
    CALL assign_if_present(y)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_integers
    INTEGER :: x(5) = 1
    INTEGER :: y(5) = (/(i, i=1, 5)/)

    CALL TAG_TEST("Test_assign_with_present_integers")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(assert_integer_array(x, y), .TRUE.)

    CALL TAG_TEST("Test_assign_no_present_integers")
    CALL assign_if_present(y)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_real64
    REAL(dp) :: x, y

    x = 1.0
    y = 2.0
    CALL TAG_TEST("Test_assign_with_present_real64")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(x, y)

    CALL TAG_TEST("Test_assign_no_present_real64")
    CALL assign_if_present(y)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_real32
    REAL(sp) :: x, y

    x = 1.0
    y = 2.0
    CALL TAG_TEST("Test_assign_with_present_real32")
    CALL assign_if_present(y, x)
    CALL ASSERT_EQUAL(x, y)

    CALL TAG_TEST("Test_assign_no_present_real32")
    CALL assign_if_present(y)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_logical_allocatable_1d
    LOGICAL :: x(5) = .TRUE.
    LOGICAL, ALLOCATABLE :: y(:), z(:)

    CALL TAG_TEST("Test_assign_with_present_logical_allocatable_1d_not_allocated")
    CALL assign_if_present_allocatable(y, x)
    CALL ASSERT_EQUAL(assert_logical_array(x, y), .TRUE.)

    ALLOCATE (z(7))
    CALL TAG_TEST("Test_assign_with_present_logical_allocatable_1d_allocated")
    CALL assign_if_present_allocatable(z, x)
    CALL ASSERT_EQUAL(assert_logical_array(x, z), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_integer_allocatable
    INTEGER :: x = 1
    INTEGER, ALLOCATABLE :: y, z

    CALL TAG_TEST("Test_assign_with_present_integer_allocatable_not_allocated")
    CALL assign_if_present_allocatable(y, x)
    CALL ASSERT_EQUAL(x, y)

    ALLOCATE (z)
    CALL TAG_TEST("Test_assign_with_present_logical_allocatable_allocated")
    CALL assign_if_present_allocatable(z, x)
    CALL ASSERT_EQUAL(x, z)
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_integer_allocatable_1d
    INTEGER :: x(5) = (/1, 2, 3, 4, 5/)
    INTEGER, ALLOCATABLE :: y(:), z(:)

    CALL TAG_TEST("Test_assign_with_present_integer_allocatable_1d_not_allocated")
    CALL assign_if_present_allocatable(y, x)
    CALL ASSERT_EQUAL(assert_integer_array(x, y), .TRUE.)

    ALLOCATE (z(7))
    CALL TAG_TEST("Test_assign_with_present_integer_allocatable_1d_allocated")
    CALL assign_if_present_allocatable(z, x)
    CALL ASSERT_EQUAL(assert_integer_array(x, z), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_real_allocatable
    REAL(dp) :: x = 1.0
    REAL(dp), ALLOCATABLE :: y, z

    CALL TAG_TEST("Test_assign_with_present_real_allocatable_not_allocated")
    CALL assign_if_present_allocatable(y, x)
    CALL ASSERT_EQUAL(x, y)

    ALLOCATE (z)
    CALL TAG_TEST("Test_assign_with_present_logical_allocatable_allocated")
    CALL assign_if_present_allocatable(z, x)
    CALL ASSERT_EQUAL(x, z)
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_real_allocatable_1d
    REAL(dp) :: x(5) = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    REAL(dp), ALLOCATABLE :: y(:), z(:)

    CALL TAG_TEST("Test_assign_with_present_real_allocatable_1d_not_allocated")
    CALL assign_if_present_allocatable(y, x)
    CALL ASSERT_EQUAL(assert_real_array(x, y), .TRUE.)

    ALLOCATE (z(7))
    CALL TAG_TEST("Test_assign_with_present_real_allocatable_1d_allocated")
    CALL assign_if_present_allocatable(z, x)
    CALL ASSERT_EQUAL(assert_real_array(x, z), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_assign_if_present_character_allocatable
    CHARACTER(len=3) :: x
    CHARACTER(len=:), ALLOCATABLE :: y

    x = 'abc'
    CALL TAG_TEST("Test_assign_if_present_character_allocatable")
    CALL assign_if_present_allocatable(y, x)
    CALL ASSERT_EQUAL(x, y)
  END SUBROUTINE

  SUBROUTINE Test_if_associated
    REAL(dp), CONTIGUOUS, POINTER :: ptr(:, :), output(:, :)
    REAL(dp), TARGET :: arr1(10, 10), arr2(10, 10)

    arr1 = 1.0
    arr2 = 2.0

    ! NAG compiler cannot determine ASSOCIATED unless the pointer is pointed
    !   to a target or if it is a NULL pointer
    ptr => NULL()
    CALL TAG_TEST("Test_if_associated_false")
    output => if_associated(ptr)
    CALL ASSERT_EQUAL(ASSOCIATED(output), .FALSE.)

    CALL TAG_TEST("Test_if_associated_false_else")
    output => if_associated(ptr, arr2)
    CALL ASSERT_EQUAL(assert_real_2d_array(output, arr2), .TRUE.)

    ptr => arr1
    CALL TAG_TEST("Test_if_associated_true")
    output => if_associated(ptr)
    CALL ASSERT_EQUAL(assert_real_2d_array(output, arr1), .TRUE.)

    CALL TAG_TEST("Test_if_associated_true_else")
    output => if_associated(ptr, arr2)
    CALL ASSERT_EQUAL(assert_real_2d_array(output, arr1), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_swap_int
    INTEGER :: a, b

    a = 1
    b = 2
    CALL TAG_TEST("Test_swap_int_a")
    CALL swap(a, b)
    CALL ASSERT_EQUAL(a, 2)
    CALL TAG_TEST("Test_swap_int_b")
    CALL ASSERT_EQUAL(b, 1)
  END SUBROUTINE

  SUBROUTINE Test_resize_arr_c1d
    CHARACTER(len=256), ALLOCATABLE :: arr(:)
    INTEGER :: nelem = 12
    CHARACTER(len=100) :: log_in_file, logfile

    CALL TAG_TEST("Test_resize_arr_c1d_not_allocated")
    CALL resize_arr_c1d(arr, nelem)
    CALL ASSERT_EQUAL(SIZE(arr), 1)

    CALL TAG_TEST("Test_resize_arr_c1d_allocated")
    CALL resize_arr_c1d(arr, nelem)
    CALL ASSERT_EQUAL(SIZE(arr), 13)

    CALL TAG_TEST("Test_resize_arr_c1d_allocated2")
    CALL resize_arr_c1d(arr, nelem)
    CALL ASSERT_EQUAL(SIZE(arr), 25)

    CALL TAG_TEST("Test_reisze_arr_c1d_nelem_error")
    logfile = 'logger_output.txt'

    CALL open_new_logfile(nerr, TRIM(logfile))
    CALL init_logger(5, .FALSE., nerr, callback_abort=custom_exit_dummy)

    CALL resize_arr_c1d(arr, -1)
    CLOSE (nerr)

    CALL open_logfile(nerr, TRIM(logfile))
    READ (nerr, '(A)') log_in_file

    CALL STRING_CONTAINS(log_in_file, "FINISH PE:     5 &
      &mo_fortran_tools:resize_arr_c1d: nelem must be > 0")
    CLOSE (nerr)
  END SUBROUTINE

  SUBROUTINE Test_copy_1d_dp
    REAL(dp) :: src(10) = 1.0, dest(10)

    CALL TAG_TEST("Test_copy_1d_dp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_1d_dp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_2d_dp
    REAL(dp) :: src(10, 10) = 1.0, dest(10, 10)

    CALL TAG_TEST("Test_copy_2d_dp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_2d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_2d_dp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_2d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_3d_dp
    REAL(dp) :: src(10, 10, 10) = 1.0, dest(10, 10, 10)

    CALL TAG_TEST("Test_copy_3d_dp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_3d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_3d_dp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_3d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_4d_dp
    REAL(dp) :: src(5, 5, 5, 5) = 1.0, dest(5, 5, 5, 5)

    CALL TAG_TEST("Test_copy_4d_dp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_4d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_4d_dp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_4d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_5d_dp
    REAL(dp) :: src(5, 5, 5, 5, 5) = 1.0, dest(5, 5, 5, 5, 5)

    CALL TAG_TEST("Test_copy_5d_dp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_5d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_5d_dp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_5d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_5d_sp
    REAL(sp) :: src(5, 5, 5, 5, 5) = 1.0, dest(5, 5, 5, 5, 5)

    CALL TAG_TEST("Test_copy_5d_sp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_5d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_5d_sp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_5d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_2d_spdp
    REAL(sp) :: src(10, 10) = 1.0
    REAL(dp) :: dest(10, 10)

    CALL TAG_TEST("Test_copy_2d_spdp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_2d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_2d_spdp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_2d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_3d_spdp
    REAL(sp) :: src(10, 10, 10) = 1.0
    REAL(dp) :: dest(10, 10, 10)

    CALL TAG_TEST("Test_copy_3d_spdp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_3d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_3d_spdp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_3d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_4d_spdp
    REAL(sp) :: src(5, 5, 5, 5) = 1.0
    REAL(dp) :: dest(5, 5, 5, 5)

    CALL TAG_TEST("Test_copy_4d_spdp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_4d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_4d_spdp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_4d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_5d_spdp
    REAL(sp) :: src(5, 5, 5, 5, 5) = 1.0
    REAL(dp) :: dest(5, 5, 5, 5, 5)

    CALL TAG_TEST("Test_copy_5d_spdp_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_5d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(src)
    CALL TAG_TEST("Test_copy_5d_spdp_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_5d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_2d_i4
    INTEGER(i4) :: src(10, 10) = 1, dest(10, 10)
    REAL(sp) :: rand(10, 10)

    CALL TAG_TEST("Test_copy_2d_i4_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_2d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(rand)
    src = 1 + FLOOR(100*rand)
    CALL TAG_TEST("Test_copy_2d_i4_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_2d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_3d_i4
    INTEGER(i4) :: src(10, 10, 10) = 1, dest(10, 10, 10)
    REAL(sp) :: rand(10, 10, 10)

    CALL TAG_TEST("Test_copy_3d_i4_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_3d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(rand)
    src = 1 + FLOOR(100*rand)
    CALL TAG_TEST("Test_copy_3d_i4_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_3d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_5d_i4
    INTEGER(i4) :: src(5, 5, 5, 5, 5) = 1, dest(5, 5, 5, 5, 5)
    REAL(sp) :: rand(5, 5, 5, 5, 5)

    CALL TAG_TEST("Test_copy_5d_i4_ones")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_5d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(rand)
    src = 1 + FLOOR(100*rand)
    CALL TAG_TEST("Test_copy_5d_i4_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_5d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_copy_5d_l
    LOGICAL :: src(5, 5, 5, 5, 5) = .TRUE., dest(5, 5, 5, 5, 5)
    REAL(sp) :: rand(5, 5, 5, 5, 5)

    CALL TAG_TEST("Test_copy_5d_l_trues")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_logical_5d_array(src, dest), .TRUE.)

    CALL RANDOM_NUMBER(rand)
    src = rand < 0.5
    CALL TAG_TEST("Test_copy_5d_l_random")
    CALL copy(src, dest, .FALSE.)
    CALL ASSERT_EQUAL(assert_logical_5d_array(src, dest), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_1d_dp
    REAL(dp) :: arr(10), zeros(10) = 0.0

    CALL TAG_TEST("Test_init_zero_1d_dp")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_1d_sp
    REAL(sp) :: arr(10), zeros(10) = 0.0

    CALL TAG_TEST("Test_init_zero_1d_sp")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_2d_dp
    REAL(dp) :: arr(10, 10), zeros(10, 10) = 0.0

    CALL TAG_TEST("Test_init_zero_2d_dp")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_2d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_2d_i4
    INTEGER(i4) :: arr(10, 10), zeros(10, 10) = 0

    CALL TAG_TEST("Test_init_zero_2d_i4")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_2d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_3d_dp
    REAL(dp) :: arr(10, 10, 10), zeros(10, 10, 10) = 0.0

    CALL TAG_TEST("Test_init_zero_3d_dp")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_3d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_3d_sp
    REAL(sp) :: arr(10, 10, 10), zeros(10, 10, 10) = 0.0

    CALL TAG_TEST("Test_init_zero_3d_sp")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_3d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_3d_i4
    INTEGER(i4) :: arr(10, 10, 10), zeros(10, 10, 10) = 0

    CALL TAG_TEST("Test_init_zero_3d_i4")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_3d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_4d_dp
    REAL(dp) :: arr(5, 5, 5, 5), zeros(5, 5, 5, 5) = 0.0

    CALL TAG_TEST("Test_init_zero_4d_dp")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_4d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_4d_sp
    REAL(sp) :: arr(5, 5, 5, 5), zeros(5, 5, 5, 5) = 0.0

    CALL TAG_TEST("Test_init_zero_4d_sp")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_4d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_4d_i4
    INTEGER(i4) :: arr(5, 5, 5, 5), zeros(5, 5, 5, 5) = 0

    CALL TAG_TEST("Test_init_zero_4d_i4")
    CALL init(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_4d_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_1d_dp
    REAL(dp) :: arr(10), ones(10) = 1.0

    CALL TAG_TEST("Test_init_1d_dp")
    CALL init(arr, 1.0_dp, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_array(arr, ones), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_2d_dp
    REAL(dp) :: arr(10, 10), ones(10, 10) = 1.0

    CALL TAG_TEST("Test_init_2d_dp")
    CALL init(arr, 1.0_dp, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_2d_array(arr, ones), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_3d_dp
    REAL(dp) :: arr(10, 10, 10), ones(10, 10, 10) = 1.0

    CALL TAG_TEST("Test_init_3d_dp")
    CALL init(arr, 1.0_dp, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_3d_array(arr, ones), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_3d_spdp
    REAL(sp) :: arr(10, 10, 10)
    REAL(dp) :: ones(10, 10, 10) = 1.0

    CALL TAG_TEST("Test_init_3d_spdp")
    CALL init(arr, 1.0_dp, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_spdp_3d_array(arr, ones), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_5d_dp
    REAL(dp) :: arr(5, 5, 5, 5, 5), ones(5, 5, 5, 5, 5) = 1.0

    CALL TAG_TEST("Test_init_5d_dp")
    CALL init(arr, 1.0_dp, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_5d_array(arr, ones), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_5d_sp
    REAL(sp) :: arr(5, 5, 5, 5, 5), ones(5, 5, 5, 5, 5) = 1.0

    CALL TAG_TEST("Test_init_5d_sp")
    CALL init(arr, 1.0, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_5d_array(arr, ones), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_5d_i4
    INTEGER(i4) :: arr(5, 5, 5, 5, 5), ones(5, 5, 5, 5, 5) = 1

    CALL TAG_TEST("Test_init_5d_i4")
    CALL init(arr, 1, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_5d_array(arr, ones), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_5d_l
    LOGICAL :: arr(5, 5, 5, 5, 5), trues(5, 5, 5, 5, 5) = .TRUE.

    CALL TAG_TEST("Test_init_5d_l")
    CALL init(arr, .TRUE., .FALSE.)
    CALL ASSERT_EQUAL(assert_logical_5d_array(arr, trues), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_var_scale_3d
    REAL(dp) :: arr(10, 10, 10) = 1.0, scale = 5.0
    REAL(dp) :: ans(10, 10, 10) = 5.0

    CALL TAG_TEST("Test_var_scale_3d")
    CALL var_scale(arr, scale, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_3d_array(arr, ans), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_var_addc_3d_dp
    REAL(dp) :: arr(10, 10, 10) = 2.0, const = 3.0
    REAL(dp) :: ans(10, 10, 10) = 5.0

    CALL TAG_TEST("Test_var_addc_3d_dp")
    CALL var_add(arr, const, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_3d_array(arr, ans), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_negative2zero_4d_dp
    REAL(dp) :: arr(5, 5, 5, 5), ans(5, 5, 5, 5)
    INTEGER :: i, j, k, l

    CALL RANDOM_NUMBER(arr)
    DO i = 1, 5
      DO j = 1, 5
        DO k = 1, 5
          DO l = 1, 5
            arr(i, j, k, l) = arr(i, j, k, l) - 0.5
            IF (arr(i, j, k, l) < 0.0) THEN
              ans(i, j, k, l) = 0.0
            ELSE
              ans(i, j, k, l) = arr(i, j, k, l)
            END IF
          END DO
        END DO
      END DO
    END DO

    CALL TAG_TEST("Test_negative2zero_4d_dp")
    CALL negative2zero(arr, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_4d_array(arr, ans), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_contiguous_dp
    REAL(dp) :: arr(10), ans(10) = 7.0

    CALL TAG_TEST("Test_init_contiguous_dp")
    CALL init_contiguous_dp(arr, 10, 7.0_dp, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_array(arr, ans), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_contiguous_sp
    REAL(sp) :: arr(10), ans(10) = 7.0

    CALL TAG_TEST("Test_init_contiguous_sp")
    CALL init_contiguous_sp(arr, 10, 7.0, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_array(arr, ans), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_contiguous_dp
    REAL(dp) :: arr(10), zeros(10) = 0.0

    CALL TAG_TEST("Test_init_zero_contiguous_dp")
    CALL init_zero_contiguous_dp(arr, 10, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_zero_contiguous_sp
    REAL(sp) :: arr(10), zeros(10) = 0.0

    CALL TAG_TEST("Test_init_zero_contiguous_sp")
    CALL init_zero_contiguous_sp(arr, 10, .FALSE.)
    CALL ASSERT_EQUAL(assert_real_sp_array(arr, zeros), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_contiguous_i4
    INTEGER(i4) :: arr(10), ans(10) = 7

    CALL TAG_TEST("Test_init_contiguous_i4")
    CALL init_contiguous_i4(arr, 10, 7, .FALSE.)
    CALL ASSERT_EQUAL(assert_integer_array(arr, ans), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_init_contiguous_l
    LOGICAL :: arr(10), ans(10) = .TRUE.

    CALL TAG_TEST("Test_init_contiguous_l")
    CALL init_contiguous_l(arr, 10, .TRUE., .FALSE.)
    CALL ASSERT_EQUAL(assert_logical_array(arr, ans), .TRUE.)
  END SUBROUTINE

  SUBROUTINE Test_minval_1d
    REAL(sp) :: rand(10)
    INTEGER :: arr(10), min

    CALL RANDOM_NUMBER(rand)
    ! Make arr an array of numbers between 1 and 100
    arr = 1 + FLOOR(100*rand)
    CALL TAG_TEST("Test_minval_1d")
    min = minval_1d(arr, .FALSE.)
    CALL ASSERT_EQUAL(min, MINVAL(arr(:)))
  END SUBROUTINE

  SUBROUTINE Test_minval_2d
    REAL(sp) :: rand(10, 10)
    INTEGER :: arr(10, 10), min

    CALL RANDOM_NUMBER(rand)
    ! Make arr an array of numbers between 1 and 100
    arr = 1 + FLOOR(100*rand)
    CALL TAG_TEST("Test_minval_2d")
    min = minval_2d(arr, .FALSE.)
    CALL ASSERT_EQUAL(min, MINVAL(arr(:, :)))
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_dp_3_2
    REAL(dp), POINTER :: ptr_out(:, :, :)
    REAL(dp), TARGET :: ptr_in(5, 10)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_dp_3_2_test2
    REAL(dp), POINTER :: ptr_out(:, :, :)
    REAL(dp), TARGET :: ptr_in(1, 10)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_dp_3_2_test3
    REAL(dp), POINTER :: ptr_out(:, :, :)
    REAL(dp), TARGET :: ptr_in(5, 1)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test3_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test3_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test3_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test3_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test3_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test3_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_dp_3_2_test4
    REAL(dp), POINTER :: ptr_out(:, :, :)
    REAL(dp), TARGET :: ptr_in(1, 1)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test4_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test4_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test4_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test4_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test4_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test4_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_dp_3_2_test5
    REAL(dp), POINTER :: ptr_out(:, :, :)
    REAL(dp), TARGET :: ptr_in(0, 0)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test5_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 0)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test5_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test5_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test5_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test5_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 0)

    CALL TAG_TEST("Test_insert_dimension_r_dp_3_2_test5_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_sp_3_2
    REAL(sp), POINTER :: ptr_out(:, :, :)
    REAL(sp), TARGET :: ptr_in(5, 10)
    ptr_in = 1.0

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_sp_3_2_test2
    REAL(sp), POINTER :: ptr_out(:, :, :)
    REAL(sp), TARGET :: ptr_in(1, 10)
    ptr_in = 1.0

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_sp_3_2_test3
    REAL(sp), POINTER :: ptr_out(:, :, :)
    REAL(sp), TARGET :: ptr_in(5, 1)
    ptr_in = 1.0

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test3_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test3_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test3_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test3_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test3_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test3_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_sp_3_2_test4
    REAL(sp), POINTER :: ptr_out(:, :, :)
    REAL(sp), TARGET :: ptr_in(1, 1)
    ptr_in = 1.0

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test4_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test4_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test4_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test4_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test4_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test4_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_sp_3_2_test5
    REAL(sp), POINTER :: ptr_out(:, :, :)
    REAL(sp), TARGET :: ptr_in(0, 0)
    ptr_in = 1.0

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test5_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 0)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test5_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test5_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test5_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test5_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 0)

    CALL TAG_TEST("Test_insert_dimension_r_sp_3_2_test5_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_i4_3_2
    INTEGER(i4), POINTER :: ptr_out(:, :, :)
    INTEGER(i4), TARGET :: ptr_in(5, 10)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_i4_3_2_test2
    INTEGER(i4), POINTER :: ptr_out(:, :, :)
    INTEGER(i4), TARGET :: ptr_in(1, 10)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_i4_3_2_test3
    INTEGER(i4), POINTER :: ptr_out(:, :, :)
    INTEGER(i4), TARGET :: ptr_in(5, 1)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test3_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test3_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test3_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test3_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test3_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test3_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_i4_3_2_test4
    INTEGER(i4), POINTER :: ptr_out(:, :, :)
    INTEGER(i4), TARGET :: ptr_in(1, 1)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test4_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test4_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test4_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test4_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test4_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test4_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_i4_3_2_test5
    INTEGER(i4), POINTER :: ptr_out(:, :, :)
    INTEGER(i4), TARGET :: ptr_in(0, 0)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test5_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 0)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test5_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test5_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test5_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test5_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 0)

    CALL TAG_TEST("Test_insert_dimension_i4_3_2_test5_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_l_3_2
    LOGICAL, POINTER :: ptr_out(:, :, :)
    LOGICAL, TARGET :: ptr_in(5, 10)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_l_3_2_test2
    LOGICAL, POINTER :: ptr_out(:, :, :)
    LOGICAL, TARGET :: ptr_in(1, 10)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test2_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test2_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test2_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test2_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test2_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test2_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 10)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_l_3_2_test3
    LOGICAL, POINTER :: ptr_out(:, :, :)
    LOGICAL, TARGET :: ptr_in(5, 1)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test3_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 5)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test3_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test3_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test3_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test3_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 5)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test3_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_l_3_2_test4
    LOGICAL, POINTER :: ptr_out(:, :, :)
    LOGICAL, TARGET :: ptr_in(1, 1)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test4_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test4_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test4_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test4_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test4_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test4_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_l_3_2_test5
    LOGICAL, POINTER :: ptr_out(:, :, :)
    LOGICAL, TARGET :: ptr_in(0, 0)

    CALL insert_dimension(ptr_out, ptr_in, 2)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test5_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 0)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test5_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test5_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)

    CALL insert_dimension(ptr_out, ptr_in, 1)
    CALL TAG_TEST("Test_insert_dimension_l_3_2_test5_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 1)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test5_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 0)

    CALL TAG_TEST("Test_insert_dimension_l_3_2_test5_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 0)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_dp_6_5
    REAL(dp), POINTER :: ptr_out(:, :, :, :, :, :)
    REAL(dp), TARGET :: ptr_in(2, 3, 4, 5, 6)

    CALL insert_dimension(ptr_out, ptr_in, 3)
    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 2)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 3)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_fourth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 4), 4)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_fifth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 5), 5)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_sixth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 6), 6)

    CALL insert_dimension(ptr_out, ptr_in, 6)
    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 2)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 3)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 4)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_fourth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 4), 5)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_fifth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 5), 6)

    CALL TAG_TEST("Test_insert_dimension_r_dp_6_5_sixth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 6), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_r_sp_6_5
    REAL(sp), POINTER :: ptr_out(:, :, :, :, :, :)
    REAL(sp), TARGET :: ptr_in(2, 3, 4, 5, 6)

    CALL insert_dimension(ptr_out, ptr_in, 3)
    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 2)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 3)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_fourth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 4), 4)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_fifth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 5), 5)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_sixth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 6), 6)

    CALL insert_dimension(ptr_out, ptr_in, 6)
    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 2)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 3)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 4)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_fourth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 4), 5)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_fifth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 5), 6)

    CALL TAG_TEST("Test_insert_dimension_r_sp_6_5_sixth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 6), 1)
  END SUBROUTINE

  SUBROUTINE Test_insert_dimension_i4_6_5
    INTEGER(i4), POINTER :: ptr_out(:, :, :, :, :, :)
    INTEGER(i4), TARGET :: ptr_in(2, 3, 4, 5, 6)

    CALL insert_dimension(ptr_out, ptr_in, 3)
    CALL TAG_TEST("Test_insert_dimension_i4_6_5_first_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 2)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_second_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 3)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_third_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 1)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_fourth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 4), 4)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_fifth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 5), 5)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_sixth_dim")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 6), 6)

    CALL insert_dimension(ptr_out, ptr_in, 6)
    CALL TAG_TEST("Test_insert_dimension_i4_6_5_first_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 1), 2)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_second_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 2), 3)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_third_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 3), 4)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_fourth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 4), 5)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_fifth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 5), 6)

    CALL TAG_TEST("Test_insert_dimension_i4_6_5_sixth_dim2")
    CALL ASSERT_EQUAL(SIZE(ptr_out, 6), 1)
  END SUBROUTINE

  SUBROUTINE Test_DO_DEALLOCATE_r4D
    REAL(dp), ALLOCATABLE :: arr(:, :, :, :)

    ALLOCATE (arr(5, 5, 5, 5))
    CALL TAG_TEST("Test_DO_DEALLOCATE_r4D")
    CALL DO_DEALLOCATE(arr)
    CALL ASSERT_EQUAL(ALLOCATED(arr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_DEALLOCATE_r3D
    REAL(dp), ALLOCATABLE :: arr(:, :, :)

    ALLOCATE (arr(10, 10, 10))
    CALL TAG_TEST("Test_DO_DEALLOCATE_r3D")
    CALL DO_DEALLOCATE(arr)
    CALL ASSERT_EQUAL(ALLOCATED(arr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_DEALLOCATE_r2D
    REAL(dp), ALLOCATABLE :: arr(:, :)

    ALLOCATE (arr(10, 10))
    CALL TAG_TEST("Test_DO_DEALLOCATE_r2D")
    CALL DO_DEALLOCATE(arr)
    CALL ASSERT_EQUAL(ALLOCATED(arr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_DEALLOCATE_r1D
    REAL(dp), ALLOCATABLE :: arr(:)

    ALLOCATE (arr(10))
    CALL TAG_TEST("Test_DO_DEALLOCATE_r1D")
    CALL DO_DEALLOCATE(arr)
    CALL ASSERT_EQUAL(ALLOCATED(arr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_DEALLOCATE_i3D
    INTEGER, ALLOCATABLE :: arr(:, :, :)

    ALLOCATE (arr(10, 10, 10))
    CALL TAG_TEST("Test_DO_DEALLOCATE_i3D")
    CALL DO_DEALLOCATE(arr)
    CALL ASSERT_EQUAL(ALLOCATED(arr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_DEALLOCATE_i2D
    INTEGER, ALLOCATABLE :: arr(:, :)

    ALLOCATE (arr(10, 10))
    CALL TAG_TEST("Test_DO_DEALLOCATE_i2D")
    CALL DO_DEALLOCATE(arr)
    CALL ASSERT_EQUAL(ALLOCATED(arr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_DEALLOCATE_i1D
    INTEGER, ALLOCATABLE :: arr(:)

    ALLOCATE (arr(10))
    CALL TAG_TEST("Test_DO_DEALLOCATE_i1D")
    CALL DO_DEALLOCATE(arr)
    CALL ASSERT_EQUAL(ALLOCATED(arr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_PTR_DEALLOCATE_r3D
    REAL(dp), POINTER :: ptr(:, :, :)

    ALLOCATE (ptr(10, 10, 10))

    CALL TAG_TEST("Test_DO_PTR_DEALLOCATE_r3D")
    CALL DO_PTR_DEALLOCATE(ptr)
    CALL ASSERT_EQUAL(ASSOCIATED(ptr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_PTR_DEALLOCATE_r2D
    REAL(dp), POINTER :: ptr(:, :)

    ALLOCATE (ptr(10, 10))

    CALL TAG_TEST("Test_DO_PTR_DEALLOCATE_r2D")
    CALL DO_PTR_DEALLOCATE(ptr)
    CALL ASSERT_EQUAL(ASSOCIATED(ptr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_PTR_DEALLOCATE_dp1D
    REAL(dp), POINTER :: ptr(:)

    ALLOCATE (ptr(10))

    CALL TAG_TEST("Test_DO_PTR_DEALLOCATE_dp1D")
    CALL DO_PTR_DEALLOCATE(ptr)
    CALL ASSERT_EQUAL(ASSOCIATED(ptr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_PTR_DEALLOCATE_sp1D
    REAL(sp), POINTER :: ptr(:)

    ALLOCATE (ptr(10))

    CALL TAG_TEST("Test_DO_PTR_DEALLOCATE_sp1D")
    CALL DO_PTR_DEALLOCATE(ptr)
    CALL ASSERT_EQUAL(ASSOCIATED(ptr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_DO_PTR_DEALLOCATE_int1D
    INTEGER, POINTER :: ptr(:)

    ALLOCATE (ptr(10))

    CALL TAG_TEST("Test_DO_PTR_DEALLOCATE_int1D")
    CALL DO_PTR_DEALLOCATE(ptr)
    CALL ASSERT_EQUAL(ASSOCIATED(ptr), .FALSE.)
  END SUBROUTINE

  SUBROUTINE Test_assert_acc_host_only
    ! OpenACC version is left TODO
    CALL TAG_TEST("Test_assert_acc_host_only_true")
    CALL assert_acc_host_only("Unit_test", .TRUE.)
    CALL SUCCEED

    CALL TAG_TEST("Test_assert_acc_host_only_false")
    CALL assert_acc_host_only("Unit_test", .FALSE.)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assert_acc_device_only
    ! OpenACC version is left TODO
    CALL TAG_TEST("Test_assert_acc_device_only_true")
    CALL assert_acc_device_only("Unit_test", .TRUE.)
    CALL SUCCEED

    CALL TAG_TEST("Test_assert_acc_device_only_false")
    CALL assert_acc_device_only("Unit_test", .FALSE.)
    CALL SUCCEED
  END SUBROUTINE

  SUBROUTINE Test_assert_lacc_equals_i_am_accel_node
    ! OpenACC version is left TODO
    CALL TAG_TEST("Test_assert_lacc_equals_i_am_accel_node_match_true")
    CALL assert_lacc_equals_i_am_accel_node("Unit_test", .TRUE., .TRUE.)
    CALL SUCCEED

    CALL TAG_TEST("Test_assert_lacc_equals_i_am_accel_node_match_false")
    CALL assert_lacc_equals_i_am_accel_node("Unit_test", .FALSE., .FALSE.)
    CALL SUCCEED

    CALL TAG_TEST("Test_assert_lacc_equals_i_am_accel_node_false")
    CALL assert_lacc_equals_i_am_accel_node("Unit_test", .FALSE., .TRUE.)
    CALL SUCCEED
  END SUBROUTINE

! Support functions for testing

  LOGICAL FUNCTION assert_logical_array(array1, array2)
    LOGICAL, INTENT(IN) :: array1(:), array2(:)
    INTEGER :: i

    assert_logical_array = .TRUE.
    DO i = 1, SIZE(array1)
      IF (array1(i) .NEQV. array2(i)) THEN
        assert_logical_array = .FALSE.
        EXIT
      END IF
    END DO
  END FUNCTION assert_logical_array

  LOGICAL FUNCTION assert_integer_array(array1, array2)
    INTEGER, INTENT(IN) :: array1(:), array2(:)
    INTEGER :: i

    assert_integer_array = .TRUE.
    DO i = 1, SIZE(array1)
      IF (array1(i) /= array2(i)) THEN
        assert_integer_array = .FALSE.
        EXIT
      END IF
    END DO
  END FUNCTION assert_integer_array

  LOGICAL FUNCTION assert_real_array(array1, array2)
    REAL(dp), INTENT(IN) :: array1(:), array2(:)
    INTEGER :: i

    assert_real_array = .TRUE.
    DO i = 1, SIZE(array1)
      IF (array1(i) /= array2(i)) THEN
        assert_real_array = .FALSE.
        EXIT
      END IF
    END DO
  END FUNCTION assert_real_array

  LOGICAL FUNCTION assert_real_2d_array(array1, array2)
    REAL(dp), INTENT(IN) :: array1(:, :), array2(:, :)
    INTEGER :: i, j

    assert_real_2d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        IF (array1(i, j) /= array2(i, j)) THEN
          assert_real_2d_array = .FALSE.
          EXIT
        END IF
      END DO
    END DO
  END FUNCTION assert_real_2d_array

  LOGICAL FUNCTION assert_real_3d_array(array1, array2)
    REAL(dp), INTENT(IN) :: array1(:, :, :), array2(:, :, :)
    INTEGER :: i, j, k

    assert_real_3d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          IF (array1(i, j, k) /= array2(i, j, k)) THEN
            assert_real_3d_array = .FALSE.
            EXIT
          END IF
        END DO
      END DO
    END DO
  END FUNCTION assert_real_3d_array

  LOGICAL FUNCTION assert_real_4d_array(array1, array2)
    REAL(dp), INTENT(IN) :: array1(:, :, :, :), array2(:, :, :, :)
    INTEGER :: i, j, k, l

    assert_real_4d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            IF (array1(i, j, k, l) /= array2(i, j, k, l)) THEN
              assert_real_4d_array = .FALSE.
              EXIT
            END IF
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_real_4d_array

  LOGICAL FUNCTION assert_real_5d_array(array1, array2)
    REAL(dp), INTENT(IN) :: array1(:, :, :, :, :), array2(:, :, :, :, :)
    INTEGER :: i, j, k, l, m

    assert_real_5d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            DO m = 1, SIZE(array1, 5)
              IF (array1(i, j, k, l, m) /= array2(i, j, k, l, m)) THEN
                assert_real_5d_array = .FALSE.
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_real_5d_array

  LOGICAL FUNCTION assert_real_sp_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:), array2(:)
    INTEGER :: i

    assert_real_sp_array = .TRUE.
    DO i = 1, SIZE(array1)
      IF (array1(i) /= array2(i)) THEN
        assert_real_sp_array = .FALSE.
        EXIT
      END IF
    END DO
  END FUNCTION assert_real_sp_array

  LOGICAL FUNCTION assert_real_sp_3d_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:, :, :), array2(:, :, :)
    INTEGER :: i, j, k

    assert_real_sp_3d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          IF (array1(i, j, k) /= array2(i, j, k)) THEN
            assert_real_sp_3d_array = .FALSE.
            EXIT
          END IF
        END DO
      END DO
    END DO
  END FUNCTION assert_real_sp_3d_array

  LOGICAL FUNCTION assert_real_sp_4d_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:, :, :, :), array2(:, :, :, :)
    INTEGER :: i, j, k, l

    assert_real_sp_4d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            IF (array1(i, j, k, l) /= array2(i, j, k, l)) THEN
              assert_real_sp_4d_array = .FALSE.
              EXIT
            END IF
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_real_sp_4d_array

  LOGICAL FUNCTION assert_real_sp_5d_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:, :, :, :, :), array2(:, :, :, :, :)
    INTEGER :: i, j, k, l, m

    assert_real_sp_5d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            DO m = 1, SIZE(array1, 5)
              IF (array1(i, j, k, l, m) /= array2(i, j, k, l, m)) THEN
                assert_real_sp_5d_array = .FALSE.
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_real_sp_5d_array

  LOGICAL FUNCTION assert_real_spdp_2d_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:, :)
    REAL(dp), INTENT(IN) :: array2(:, :)
    INTEGER :: i, j

    assert_real_spdp_2d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        IF (array1(i, j) /= array2(i, j)) THEN
          assert_real_spdp_2d_array = .FALSE.
          EXIT
        END IF
      END DO
    END DO
  END FUNCTION assert_real_spdp_2d_array

  LOGICAL FUNCTION assert_real_spdp_3d_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:, :, :)
    REAL(dp), INTENT(IN) :: array2(:, :, :)
    INTEGER :: i, j, k

    assert_real_spdp_3d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          IF (array1(i, j, k) /= array2(i, j, k)) THEN
            assert_real_spdp_3d_array = .FALSE.
            EXIT
          END IF
        END DO
      END DO
    END DO
  END FUNCTION assert_real_spdp_3d_array

  LOGICAL FUNCTION assert_real_spdp_4d_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:, :, :, :)
    REAL(dp), INTENT(IN) :: array2(:, :, :, :)
    INTEGER :: i, j, k, l

    assert_real_spdp_4d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            IF (array1(i, j, k, l) /= array2(i, j, k, l)) THEN
              assert_real_spdp_4d_array = .FALSE.
              EXIT
            END IF
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_real_spdp_4d_array

  LOGICAL FUNCTION assert_real_spdp_5d_array(array1, array2)
    REAL(sp), INTENT(IN) :: array1(:, :, :, :, :)
    REAL(dp), INTENT(IN) :: array2(:, :, :, :, :)
    INTEGER :: i, j, k, l, m

    assert_real_spdp_5d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            DO m = 1, SIZE(array1, 5)
              IF (array1(i, j, k, l, m) /= array2(i, j, k, l, m)) THEN
                assert_real_spdp_5d_array = .FALSE.
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_real_spdp_5d_array

  LOGICAL FUNCTION assert_integer_2d_array(array1, array2)
    INTEGER(i4), INTENT(IN) :: array1(:, :), array2(:, :)
    INTEGER :: i, j

    assert_integer_2d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        IF (array1(i, j) /= array2(i, j)) THEN
          assert_integer_2d_array = .FALSE.
          EXIT
        END IF
      END DO
    END DO
  END FUNCTION assert_integer_2d_array

  LOGICAL FUNCTION assert_integer_3d_array(array1, array2)
    INTEGER(i4), INTENT(IN) :: array1(:, :, :), array2(:, :, :)
    INTEGER :: i, j, k

    assert_integer_3d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          IF (array1(i, j, k) /= array2(i, j, k)) THEN
            assert_integer_3d_array = .FALSE.
            EXIT
          END IF
        END DO
      END DO
    END DO
  END FUNCTION assert_integer_3d_array

  LOGICAL FUNCTION assert_integer_4d_array(array1, array2)
    INTEGER(i4), INTENT(IN) :: array1(:, :, :, :), array2(:, :, :, :)
    INTEGER :: i, j, k, l

    assert_integer_4d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            IF (array1(i, j, k, l) /= array2(i, j, k, l)) THEN
              assert_integer_4d_array = .FALSE.
              EXIT
            END IF
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_integer_4d_array

  LOGICAL FUNCTION assert_integer_5d_array(array1, array2)
    INTEGER(i4), INTENT(IN) :: array1(:, :, :, :, :), array2(:, :, :, :, :)
    INTEGER :: i, j, k, l, m

    assert_integer_5d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            DO m = 1, SIZE(array1, 5)
              IF (array1(i, j, k, l, m) /= array2(i, j, k, l, m)) THEN
                assert_integer_5d_array = .FALSE.
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_integer_5d_array

  LOGICAL FUNCTION assert_logical_5d_array(array1, array2)
    LOGICAL, INTENT(IN) :: array1(:, :, :, :, :), array2(:, :, :, :, :)
    INTEGER :: i, j, k, l, m

    assert_logical_5d_array = .TRUE.
    DO i = 1, SIZE(array1, 1)
      DO j = 1, SIZE(array1, 2)
        DO k = 1, SIZE(array1, 3)
          DO l = 1, SIZE(array1, 4)
            DO m = 1, SIZE(array1, 5)
              IF (array1(i, j, k, l, m) .NEQV. array2(i, j, k, l, m)) THEN
                assert_logical_5d_array = .FALSE.
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
  END FUNCTION assert_logical_5d_array

END MODULE
