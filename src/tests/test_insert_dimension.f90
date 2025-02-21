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

PROGRAM test_insert_dimension
  USE mo_kind, ONLY: wp
  USE mo_fortran_tools, ONLY: insert_dimension
  USE iso_c_binding, ONLY: c_double
  USE mo_util_system, ONLY: util_exit
  USE mo_util_timer,  ONLY: util_gettimeofday
  IMPLICIT NONE
  INTEGER, PARAMETER :: num_test_iterations=1000
  INTEGER :: test_iteration
  CALL init_rng
  DO test_iteration = 1, num_test_iterations
    CALL test_wp
    CALL test_i
  END DO
CONTAINS
  SUBROUTINE test_wp
    REAL(wp), TARGET, ALLOCATABLE :: a_wp(:, :)
    REAL(wp), POINTER :: p_wp(:, :, :)
    INTEGER :: a_shape_base(2), a_shape_use(2), a_offset(2), &
         a_stride(2), a_last_idx(2), dim_insert_pos
    INTEGER :: p_shape(3), slice_shape(2)
    REAL :: rand_val(2, 5)
    INTEGER, PARAMETER :: max_shape(2) = (/ 100, 1000 /)
    CALL RANDOM_NUMBER(rand_val)
    a_shape_base = max_shape * rand_val(:,1) + 1
    ALLOCATE(a_wp(a_shape_base(1), a_shape_base(2)))
    a_shape_use = a_shape_base * (0.5* rand_val(:,2) + 0.5) + 1
    a_offset = (a_shape_base - a_shape_use) * rand_val(:,3) + 1
    dim_insert_pos = rand_val(1,4) * 3 + 1
    IF (rand_val(2,4) >= 0.5) THEN
      a_stride = a_shape_use * rand_val(:,5) + 1
    ELSE
      a_stride = 1
    END IF
    a_last_idx = a_offset + a_shape_use - 1
    CALL RANDOM_NUMBER(a_wp)
    CALL insert_dimension(p_wp, a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
         &                           a_offset(2):a_last_idx(2):a_stride(2)), &
         &                dim_insert_pos)
    p_shape = SHAPE(p_wp)
    slice_shape = SHAPE(a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
         &                   a_offset(2):a_last_idx(2):a_stride(2)))
    SELECT CASE (dim_insert_pos)
    CASE (1)
      IF (ANY(SHAPE(p_wp(1, :, :)) /= slice_shape)) &
           CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_wp(1, :, :) /= a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
           &                        a_offset(2):a_last_idx(2):a_stride(2)))) &
           CALL incorrect_extraction
    CASE (2)
      IF (ANY(SHAPE(p_wp(:, 1, :)) /= slice_shape)) &
           CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_wp(:, 1, :) /= a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
           &                        a_offset(2):a_last_idx(2):a_stride(2)))) &
           CALL incorrect_extraction
    CASE (3)
      IF (ANY(SHAPE(p_wp(:, :, 1))  /= slice_shape)) &
           CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_wp(:, :, 1) /= a_wp(a_offset(1):a_last_idx(1):a_stride(1), &
           &                        a_offset(2):a_last_idx(2):a_stride(2)))) &
           CALL incorrect_extraction
    END SELECT
  END SUBROUTINE test_wp

  SUBROUTINE test_i
    REAL, TARGET, ALLOCATABLE :: a_r(:, :)
    INTEGER, TARGET, ALLOCATABLE :: a_i(:,:)
    INTEGER, POINTER :: p_i(:, :, :)
    INTEGER :: a_shape_base(2), a_shape_use(2), a_offset(2), &
         a_stride(2), a_last_idx(2), dim_insert_pos
    INTEGER :: p_shape(3), slice_shape(2)
    REAL :: rand_val(2, 5)
    INTEGER, PARAMETER :: max_shape(2) = (/ 200, 1000 /)
    CALL RANDOM_NUMBER(rand_val)
    a_shape_base = max_shape * rand_val(:,1) + 1
    ALLOCATE(a_i(a_shape_base(1), a_shape_base(2)), &
         &   a_r(a_shape_base(1), a_shape_base(2)))
    a_shape_use = a_shape_base * (0.5* rand_val(:,2) + 0.5) + 1
    a_offset = (a_shape_base - a_shape_use) * rand_val(:,3) + 1
    dim_insert_pos = rand_val(1,4) * 3 + 1
    IF (rand_val(2,4) >= 0.5) THEN
      a_stride = a_shape_use * rand_val(:,5) + 1
    ELSE
      a_stride = 1
    END IF
    a_last_idx = a_offset + a_shape_use - 1
    CALL RANDOM_NUMBER(a_r)
    a_i = NINT(((a_r * 2.0) - 1.0) * HUGE(a_i))
    CALL insert_dimension(p_i, a_i(a_offset(1):a_last_idx(1):a_stride(1), &
         &                         a_offset(2):a_last_idx(2):a_stride(2)), &
         &                dim_insert_pos)
    p_shape = SHAPE(p_i)
    slice_shape = SHAPE(a_i(a_offset(1):a_last_idx(1):a_stride(1), &
         &                  a_offset(2):a_last_idx(2):a_stride(2)))
    SELECT CASE (dim_insert_pos)
    CASE (1)
      IF (ANY(SHAPE(p_i(1, :, :)) /= slice_shape)) &
           CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_i(1, :, :) /= a_i(a_offset(1):a_last_idx(1):a_stride(1), &
           &                      a_offset(2):a_last_idx(2):a_stride(2)))) &
           CALL incorrect_extraction
    CASE (2)
      IF (ANY(SHAPE(p_i(:, 1, :)) /= slice_shape)) &
           CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_i(:, 1, :) /= a_i(a_offset(1):a_last_idx(1):a_stride(1), &
           &                      a_offset(2):a_last_idx(2):a_stride(2)))) &
           CALL incorrect_extraction
    CASE (3)
      IF (ANY(SHAPE(p_i(:, :, 1))  /= slice_shape)) &
           CALL shape_problem(p_shape, slice_shape, dim_insert_pos)
      IF (ANY(p_i(:, :, 1) /= a_i(a_offset(1):a_last_idx(1):a_stride(1), &
           &                      a_offset(2):a_last_idx(2):a_stride(2)))) &
           CALL incorrect_extraction
    END SELECT
  END SUBROUTINE test_i

  SUBROUTINE init_rng
    INTEGER :: rseed_size
    INTEGER, ALLOCATABLE :: rseed(:)
    REAL(c_double) :: unix_time
    CALL RANDOM_SEED(size=rseed_size)
    ALLOCATE(rseed(rseed_size))
    rseed = 4711
    unix_time = util_gettimeofday()
    rseed(1) = IEOR(INT(unix_time), &
         &          INT((unix_time - FLOOR(unix_time)) * 1000000.0))
    CALL RANDOM_SEED(put=rseed)
    WRITE (0, '(a,i0)') 'rseed=', rseed(1)
  END SUBROUTINE init_rng

  SUBROUTINE shape_problem(p_shape, slice_shape, dim_insert_pos)
    INTEGER, INTENT(in) :: p_shape(3), slice_shape(2), dim_insert_pos
    WRITE(0, '(a,i0)') 'dim_insert_pos=', dim_insert_pos
    WRITE(0, '(a,"(",i0,", ", i0, ", ", i0, ")")') 'p_shape=', p_shape
    WRITE(0, '(a,"(",i0,", ", i0, ")")') 'slice_shape=', slice_shape
    WRITE(0, '(a)') 'incorrect shape'
    FLUSH(0)
    CALL util_exit(1)
  END SUBROUTINE shape_problem
  SUBROUTINE incorrect_extraction
    WRITE(0, '(a)') 'incorrect extraction'
    FLUSH(0)
    CALL util_exit(1)
  END SUBROUTINE incorrect_extraction
END PROGRAM test_insert_dimension
