! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

#define NOP(x) associate( x => x ); end associate

PROGRAM main

  USE utest
  USE yac
  USE mpi

  IMPLICIT NONE

  INTEGER :: global_rank, global_size
  LOGICAL :: is_target, is_empty_source

  CHARACTER (LEN=*), PARAMETER :: src_comp_name = 'source_comp'
  CHARACTER (LEN=*), PARAMETER :: tgt_comp_name = 'target_comp'
  CHARACTER (LEN=*), PARAMETER :: src_grid_name = 'source_grid'
  CHARACTER (LEN=*), PARAMETER :: tgt_grid_name = 'target_grid'
  CHARACTER (LEN=*), PARAMETER :: even_mask_name = 'even_mask'
  CHARACTER (LEN=*), PARAMETER :: odd_mask_name = 'odd_mask'
  TYPE(yac_string) :: mask_type_str(3)
  TYPE(yac_string) :: with_str(2)
  TYPE(yac_string) :: from_file_str(2)

  INTEGER, PARAMETER :: COLLECTION_SIZE = 3
  INTEGER, PARAMETER :: NUM_POINTSETS = 1
  INTEGER, PARAMETER :: NUM_POINTS = 9

  INTEGER :: config_from_file, with_field_mask
  INTEGER :: source_mask_type, target_mask_type

  INTEGER :: instance_id
  INTEGER :: comp_id, grid_id, point_id, field_ids(3,3,2,2)
  INTEGER :: default_mask_id, dummy_mask_id

  INTEGER :: interp_stack_config
  CHARACTER(LEN=YAC_MAX_CHARLEN) :: field_name
  TYPE(yac_string) :: src_mask_names(1)
  CHARACTER(LEN=YAC_MAX_CHARLEN) :: tgt_mask_name

  DOUBLE PRECISION :: send_field(NUM_POINTS, NUM_POINTSETS, COLLECTION_SIZE)
  DOUBLE PRECISION, TARGET :: recv_field(NUM_POINTS, COLLECTION_SIZE)
  TYPE(yac_dble_ptr) :: recv_field_ptr(COLLECTION_SIZE)

  LOGICAL :: check_recv_field

  INTEGER :: i, j, k, t, info, ierror

  INTEGER, PARAMETER :: max_opt_arg_len = 1024
  CHARACTER(max_opt_arg_len) :: config_dir
  INTEGER :: arg_len
  TYPE(yac_dble_ptr) :: buffer_ptr(NUM_POINTSETS, COLLECTION_SIZE)
  DOUBLE PRECISION, ALLOCATABLE, TARGET :: buffer(:)
  ALLOCATE(buffer(0))

  mask_type_str(1)%string = "even"
  mask_type_str(2)%string = "odd"
  mask_type_str(3)%string = "none"
  with_str(1)%string = "without"
  with_str(2)%string = "with"
  from_file_str(1)%string = "manual"
  from_file_str(2)%string = "yaml"

  ! ===================================================================

  CALL start_test('dummy_coupling9')

  CALL MPI_Init(ierror)

  CALL test(COMMAND_ARGUMENT_COUNT() == 1)
  CALL GET_COMMAND_ARGUMENT(1, config_dir, arg_len)

  CALL MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierror)
  CALL MPI_Comm_size(MPI_COMM_WORLD, global_size, ierror)

  IF (global_size /= 3) THEN
    WRITE ( * , * ) "Wrong number of processes (should be 3)"
    CALL error_exit
  ENDIF

  is_target = global_rank == 1
  is_empty_source = global_rank == 2

  IF (is_target) THEN
    CALL yac_finit ()
  ELSE
    CALL yac_finit(instance_id)
  END IF
  CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
  IF (is_target) THEN
    CALL yac_fread_config_yaml(TRIM(config_dir) // "coupling_test9.yaml")
  ELSE
    CALL yac_fread_config_yaml( &
      instance_id, TRIM(config_dir) // "coupling_test9.yaml")
  END IF

  ! define local component
  IF (is_target) THEN
    CALL yac_fdef_comp( tgt_comp_name, comp_id)
  ELSE
    CALL yac_fdef_comp( instance_id, src_comp_name, comp_id)
  END IF

  ! define grid (both components use an identical grids)
  IF (is_empty_source) THEN
     ! The empty source process defines an empty grid
     CALL yac_fdef_grid(src_grid_name, &
          0, 0, 0, [INTEGER :: ], [REAL :: ], &
          [REAL :: ], [INTEGER :: ], grid_id)
     ! define points at the vertices of the grid
     CALL yac_fdef_points(                    &
          grid_id, 0, YAC_LOCATION_CORNER, &
          [REAL :: ], [REAL :: ], point_id)

     ! define masks for vertices
     CALL yac_fdef_mask(                    &
          grid_id, 0, YAC_LOCATION_CORNER, &
          [INTEGER :: ], default_mask_id)
     CALL yac_fdef_mask_named(              &
          grid_id, 0, YAC_LOCATION_CORNER, &
          [INTEGER :: ], even_mask_name, dummy_mask_id)
     CALL yac_fdef_mask_named(              &
          grid_id, 0, YAC_LOCATION_CORNER, &
          [LOGICAL :: ], odd_mask_name, dummy_mask_id)
  ELSE
     CALL yac_fdef_grid(                               &
          MERGE(tgt_grid_name, src_grid_name, is_target), &
          (/3,3/), (/0,0/), (/0.0,1.0,2.0/), (/0.0,1.0,2.0/), grid_id)
     ! define points at the vertices of the grid
     CALL yac_fdef_points(                    &
          grid_id, (/3,3/), YAC_LOCATION_CORNER, &
          (/0.0,1.0,2.0/), (/0.0,1.0,2.0/), point_id)

     ! define masks for vertices
     CALL yac_fdef_mask(                    &
          grid_id, 3 * 3, YAC_LOCATION_CORNER, &
          (/1,0,0, 0,0,0, 0,0,0/), default_mask_id)
     CALL yac_fdef_mask_named(              &
          grid_id, 3 * 3, YAC_LOCATION_CORNER, &
          (/0,1,0, 1,0,1, 0,1,0/), even_mask_name, dummy_mask_id)
     CALL yac_fdef_mask_named(              &
          grid_id, 3 * 3, YAC_LOCATION_CORNER, &
          (/.TRUE.,.FALSE.,.TRUE., .FALSE.,.TRUE.,.FALSE., .TRUE.,.FALSE.,.TRUE./), &
          odd_mask_name, dummy_mask_id)
  ENDIF

  ! define field
  DO config_from_file = 1, 2
    DO with_field_mask = 1, 2
      DO source_mask_type = 1, 3
        DO target_mask_type = 1, 3
          WRITE (field_name, "('src2tgt_',A,'_',A,'_field_mask_',A,'_src_mask_',A,'_tgt_mask')") &
            TRIM(from_file_str(config_from_file)%string), &
            TRIM(with_str(with_field_mask)%string), &
            TRIM(mask_type_str(source_mask_type)%string), &
            TRIM(mask_type_str(target_mask_type)%string)
          IF (with_field_mask == 1) THEN
            CALL yac_fdef_field(                                &
              field_name, comp_id, (/point_id/), NUM_POINTSETS, &
              COLLECTION_SIZE, "1", YAC_TIME_UNIT_SECOND,       &
              field_ids(target_mask_type, source_mask_type,     &
                         with_field_mask, config_from_file))
          ELSE
            CALL yac_fdef_field_mask(                                    &
              field_name, comp_id, (/point_id/), (/default_mask_id/),    &
              NUM_POINTSETS, COLLECTION_SIZE, "1", YAC_TIME_UNIT_SECOND, &
              field_ids(target_mask_type, source_mask_type,              &
                        with_field_mask, config_from_file))
          END IF
        END DO
      END DO
    END DO
  END DO

  ! define manual couples
  CALL yac_fget_interp_stack_config(interp_stack_config)
  CALL yac_fadd_interp_stack_config_nnn( &
    interp_stack_config, YAC_NNN_AVG, 1, 0.0D0, 0.0D0)
  DO with_field_mask = 1, 2
    DO source_mask_type = 1, 3
      IF (source_mask_type == 1) THEN
        src_mask_names(1)%string = even_mask_name
      ELSE IF (source_mask_type == 2) THEN
        src_mask_names(1)%string = odd_mask_name
      END IF
      DO target_mask_type = 1, 3
        IF (target_mask_type == 1) THEN
          tgt_mask_name = even_mask_name
        ELSE IF (target_mask_type == 2) THEN
          tgt_mask_name = odd_mask_name
        END IF
        WRITE (field_name, "('src2tgt_manual_',A,'_field_mask_',A,'_src_mask_',A,'_tgt_mask')") &
          TRIM(with_str(with_field_mask)%string), &
          TRIM(mask_type_str(source_mask_type)%string), &
          TRIM(mask_type_str(target_mask_type)%string)
        IF (source_mask_type /= 3) THEN
          IF (target_mask_type /= 3) THEN
            IF (is_target) THEN
              CALL yac_fdef_couple (                                &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0,                          &
                src_mask_names = src_mask_names,                    &
                tgt_mask_name = tgt_mask_name )
            ELSE
              CALL yac_fdef_couple (instance_id,                    &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0,                          &
                src_mask_names = src_mask_names,                    &
                tgt_mask_name = tgt_mask_name )
            END IF
          ELSE
            IF (is_target) THEN
              CALL yac_fdef_couple (                                &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0,                          &
                src_mask_names = src_mask_names)
            ELSE
              CALL yac_fdef_couple (instance_id,                    &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0,                          &
                src_mask_names = src_mask_names)
            END IF
          END IF
        ELSE
          IF (target_mask_type /= 3) THEN
            IF (is_target) THEN
              CALL yac_fdef_couple (                                &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0,                          &
                tgt_mask_name = tgt_mask_name )
            ELSE
              CALL yac_fdef_couple (instance_id,                    &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0,                          &
                tgt_mask_name = tgt_mask_name )
            END IF
          ELSE
            IF (is_target) THEN
              CALL yac_fdef_couple (                                &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0)
            ELSE
              CALL yac_fdef_couple (instance_id,                    &
                src_comp_name, src_grid_name, field_name,           &
                tgt_comp_name, tgt_grid_name, field_name,           &
                '1', YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, &
                interp_stack_config, 0, 0)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  CALL yac_ffree_interp_stack_config(interp_stack_config)

  IF (is_target) THEN
    CALL yac_fenddef()
  ELSE
    CALL yac_fenddef(instance_id)
  END IF

  send_field =                        &
    RESHAPE(                          &
      (/ 1, 2, 3, 4, 5, 6, 7, 8, 9,   &
        11,12,13,14,15,16,17,18,19,   &
        21,22,23,24,25,26,27,28,29/), &
      (/NUM_POINTS, NUM_POINTSETS, COLLECTION_SIZE/))

  DO j = 1, COLLECTION_SIZE
    recv_field_ptr(j)%p => recv_field(:,j)
  END DO

  ! initalize buffer_ptr with empty pointers (only used for empty_source)
  DO i=1,NUM_POINTSETS
     DO j=1,COLLECTION_SIZE
        buffer_ptr(i,j)%p => buffer
     END DO
  END DO

  ! do time steps
  DO t = 1, 4

    IF (.NOT. is_target) THEN

      DO config_from_file = 1, 2
        DO with_field_mask = 1, 2
          DO source_mask_type = 1, 3
            DO target_mask_type = 1, 3
              IF(is_empty_source)THEN
                CALL yac_fput(                                      &
                      field_ids(target_mask_type, source_mask_type, &
                      with_field_mask, config_from_file),           &
                      NUM_POINTSETS, COLLECTION_SIZE,               &
                      buffer_ptr, info, ierror)
              ELSE
                CALL yac_fput(                                      &
                      field_ids(target_mask_type, source_mask_type, &
                      with_field_mask, config_from_file),           &
                      NUM_POINTS, NUM_POINTSETS, COLLECTION_SIZE,   &
                      send_field, info, ierror)
              ENDIF
              CALL yac_fwait(                                      &
                     field_ids(target_mask_type, source_mask_type, &
                               with_field_mask, config_from_file))
              CALL test(ierror == 0)
              CALL test(info == YAC_ACTION_COUPLING)

            END DO
          END DO
        END DO
      END DO

    ELSE

      DO config_from_file = 1, 2
        DO with_field_mask = 1, 2
          DO source_mask_type = 1, 3
            DO target_mask_type = 1, 3

              ! initialise recv_field
              recv_field = -1

              IF (IAND(t, 1) == 1) THEN
                CALL yac_fget(                                  &
                  field_ids(target_mask_type, source_mask_type, &
                            with_field_mask, config_from_file), &
                  NUM_POINTS, COLLECTION_SIZE,                  &
                  recv_field, info, ierror)
              ELSE
                CALL yac_fget_async(                            &
                  field_ids(target_mask_type, source_mask_type, &
                            with_field_mask, config_from_file), &
                  COLLECTION_SIZE, recv_field_ptr, info, ierror)
                CALL yac_fwait(                                 &
                  field_ids(target_mask_type, source_mask_type, &
                            with_field_mask, config_from_file))
              END IF
              CALL test(ierror == 0)
              CALL test(info == YAC_ACTION_COUPLING)

              DO j = 1, COLLECTION_SIZE
                DO k = 1, NUM_POINTS
                  check_recv_field =                               &
                    check_tgt(                                     &
                      recv_field(k,j), k, j, with_field_mask == 2, &
                      target_mask_type, source_mask_type)
                  CALL test(.NOT. check_recv_field)
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO

    END IF

  END DO

  IF (is_target) THEN
    CALL yac_ffinalize()
  ELSE
    CALL yac_ffinalize(instance_id)
  END IF

  CALL MPI_Finalize(ierror)

  CALL stop_test
  CALL exit_tests

! ----------------------------------------------------------

CONTAINS

  SUBROUTINE error_exit ()

    USE mpi, ONLY : mpi_abort, MPI_COMM_WORLD
    USE utest

    INTEGER :: ierror

    CALL test ( .FALSE. )
    CALL stop_test
    CALL exit_tests
    CALL mpi_abort ( MPI_COMM_WORLD, 999, ierror )

  END SUBROUTINE error_exit

  FUNCTION check_src_even_mask(value, idx, collection_idx, with_field_mask)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: collection_idx
    LOGICAL, INTENT(IN) :: with_field_mask
    LOGICAL :: check_src_even_mask

    NOP(idx)
    NOP(collection_idx)
    NOP(with_field_mask)
    check_src_even_mask = (value < 0.0) .OR. (IAND(INT(value), 1) /= 0)
  END FUNCTION

  FUNCTION check_src_odd_mask(value, idx, collection_idx, with_field_mask)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: collection_idx
    LOGICAL, INTENT(IN) :: with_field_mask
    LOGICAL :: check_src_odd_mask

    NOP(idx)
    NOP(collection_idx)
    NOP(with_field_mask)
    check_src_odd_mask =  (value < 0.0) .OR. (IAND(INT(value), 1) /= 1)
  END FUNCTION

  FUNCTION check_src_no_mask(value, idx, collection_idx, with_field_mask)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: collection_idx
    LOGICAL, INTENT(IN) :: with_field_mask
    LOGICAL :: check_src_no_mask

    check_src_no_mask = &
      MERGE( &
        INT(value) /= 1 + 10 * (collection_idx - 1), &
        INT(value) /= idx + 10 * (collection_idx - 1), with_field_mask)
  END FUNCTION

  FUNCTION check_tgt_even_mask( &
    value, idx, with_field_mask, check_src)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    LOGICAL, INTENT(IN) :: with_field_mask
    LOGICAL :: check_src
    LOGICAL :: check_tgt_even_mask

    NOP(with_field_mask)
    check_tgt_even_mask = MERGE(check_src, value /= -1.0d0, IAND(idx, 1) == 0)
  END FUNCTION

  FUNCTION check_tgt_odd_mask( &
    value, idx, with_field_mask, check_src)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    LOGICAL, INTENT(IN) :: with_field_mask
    LOGICAL :: check_src
    LOGICAL :: check_tgt_odd_mask

    NOP(with_field_mask)
    check_tgt_odd_mask = MERGE(check_src, value /= -1.0d0, IAND(idx, 1) == 1)
  END FUNCTION

  FUNCTION check_tgt_no_mask( &
    value, idx, with_field_mask, check_src)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    LOGICAL, INTENT(IN) :: with_field_mask
    LOGICAL :: check_src
    LOGICAL :: check_tgt_no_mask

    check_tgt_no_mask = &
      MERGE(value /= -1.0d0, check_src, with_field_mask .AND. (idx /= 1))
  END FUNCTION

  FUNCTION check_src( &
    value, idx, collection_idx, with_field_mask, source_mask_type)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: collection_idx
    LOGICAL, INTENT(IN) :: with_field_mask
    INTEGER :: source_mask_type
    LOGICAL :: check_src

    SELECT CASE (source_mask_type)
      CASE (1)
        check_src =                                        &
          check_src_even_mask(                             &
            value, idx, collection_idx, with_field_mask)
      CASE (2)
        check_src =                                        &
          check_src_odd_mask(                              &
            value, idx, collection_idx, with_field_mask)
      CASE (3)
        check_src =                                        &
          check_src_no_mask(                               &
            value, idx, collection_idx, with_field_mask)
      CASE DEFAULT
        check_src = .TRUE.
    END SELECT
  END FUNCTION

  FUNCTION check_tgt( &
    value, idx, collection_idx, with_field_mask, &
    target_mask_type, source_mask_type)
    DOUBLE PRECISION, INTENT(IN) :: value
    INTEGER, INTENT(IN) :: idx
    INTEGER, INTENT(IN) :: collection_idx
    LOGICAL, INTENT(IN) :: with_field_mask
    INTEGER :: target_mask_type
    INTEGER :: source_mask_type
    LOGICAL :: check_tgt

    SELECT CASE (target_mask_type)
      CASE (1)
        check_tgt =                                        &
          check_tgt_even_mask(                             &
            value, idx, with_field_mask,                   &
            check_src(                                     &
              value, idx, collection_idx, with_field_mask, &
              source_mask_type))
      CASE (2)
        check_tgt =                                        &
          check_tgt_odd_mask(                              &
            value, idx, with_field_mask,                   &
            check_src(                                     &
              value, idx, collection_idx, with_field_mask, &
              source_mask_type))
      CASE (3)
        check_tgt =                                        &
          check_tgt_no_mask(                               &
            value, idx, with_field_mask,                   &
            check_src(                                     &
              value, idx, collection_idx, with_field_mask, &
              source_mask_type))
      CASE DEFAULT
        check_tgt = .TRUE.
    END SELECT
  END FUNCTION

END PROGRAM main
