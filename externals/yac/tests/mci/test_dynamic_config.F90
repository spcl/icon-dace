! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

#define YAC_RAD (0.01745329251994329576923690768489d0) ! M_PI / 180

PROGRAM main

  use mpi
  use utest
  use yac
  use, intrinsic :: iso_c_binding, only : c_null_char

  IMPLICIT NONE

  INTEGER :: rank, global_size, ierror
  INTEGER, PARAMETER :: n_cells(2) = (/360, 180/)
  INTEGER, PARAMETER :: n_corners(2) = (/361, 181/)
  INTEGER, PARAMETER :: cyclic(2) = (/0, 0/)
  DOUBLE PRECISION, PARAMETER :: delta_lon=1., delta_lat=1.
  DOUBLE PRECISION, ALLOCATABLE :: data(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: cell_lon(:), cell_lat(:), corner_lon(:), corner_lat(:)
  INTEGER :: i, comp_id, grid_id, points_id, info, timestep_counter
  INTEGER :: field_1, field_2, field_3, field_4
  CHARACTER(len=YAC_MAX_CHARLEN), PARAMETER :: &
    weight_file_name = "test_dynamic_config.nc"
  CHARACTER(len=YAC_MAX_CHARLEN), PARAMETER :: &
    weight_file_name_dummy = "test_dynamic_config_dummy.nc"
  CHARACTER(len=YAC_MAX_CHARLEN), PARAMETER :: &
    test_config_filename_yaml = "test_dynamic_config.yaml"
  CHARACTER(len=YAC_MAX_CHARLEN), PARAMETER :: &
    test_config_filename_json = "test_dynamic_config.json"
  CHARACTER(len=YAC_MAX_CHARLEN), PARAMETER :: &
    grid_output_filename = "test_dynamic_config_grids.nc"
  LOGICAL :: file_exists

  INTERFACE
    SUBROUTINE C_UNLINK ( path ) BIND ( c, name='unlink' )
      USE, INTRINSIC :: iso_c_binding, only : c_char
      CHARACTER(KIND=c_char), DIMENSION(*) :: path
    END SUBROUTINE C_UNLINK
  END INTERFACE

  CALL MPI_Init(ierror)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)
  CALL MPI_Comm_size(MPI_COMM_WORLD, global_size, ierror)

  IF (global_size /= 4) THEN
    WRITE ( * , * ) "Wrong number of processes (should be 4)"
    CALL error_exit
  ENDIF

  ALLOCATE(data(n_cells(1), n_cells(2)))
  ALLOCATE(corner_lon(n_corners(1)), corner_lat(n_corners(2)))
  ALLOCATE(cell_lon(n_cells(1)), cell_lat(n_cells(2)))

  DO i = 1, n_corners(1)
     corner_lon(i) = (DBLE(i-1)*delta_lon-180.)*YAC_RAD
  END DO
  DO i = 1, n_corners(2)
     corner_lat(i) = (DBLE(i-1)*delta_lat-90.)*YAC_RAD
  END DO
  DO i = 1, n_cells(1)
     cell_lon(i) = ((DBLE(i)-0.5)*delta_lon-180.)*YAC_RAD
  END DO
  DO i = 1, n_cells(2)
     cell_lat(i) = ((DBLE(i)-0.5)*delta_lat-90.)*YAC_RAD
  END DO

  IF (rank == 0) CALL gen_weight_file()

  timestep_counter = 0

  IF ( MODULO(rank, 2) .eq. 0 ) THEN
     CALL comp1_main()
  ELSE
     CALL comp2_main()
  END IF

  IF (rank == 0) THEN
    INQUIRE(FILE=test_config_filename_yaml, EXIST=file_exists)
    CALL test(file_exists)
    INQUIRE(FILE=test_config_filename_json, EXIST=file_exists)
    CALL test(file_exists)
    INQUIRE(FILE=grid_output_filename, EXIST=file_exists)
    CALL test(file_exists)
    CALL C_UNLINK(TRIM(grid_output_filename) // c_null_char)
    CALL C_UNLINK(TRIM(test_config_filename_json) // c_null_char)
    CALL C_UNLINK(TRIM(test_config_filename_yaml) // c_null_char)
    CALL C_UNLINK(TRIM(weight_file_name) // c_null_char)
    CALL C_UNLINK(TRIM(weight_file_name_dummy) // c_null_char)
  END IF

  DEALLOCATE(data, cell_lon, cell_lat, corner_lon, corner_lat)

  IF (timestep_counter .NE. 11) STOP 2

  CALL MPI_Finalize(ierror)

CONTAINS

  SUBROUTINE gen_weight_file()

    INTEGER :: yac_id
    INTEGER :: grid_id_out, grid_id_in
    INTEGER :: point_id_in, point_id_out
    INTEGER :: field_id_in, field_id_out
    INTEGER :: interp_stack

    CALL yac_finit_comm(MPI_COMM_SELF, yac_id)
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fdef_datetime( &
      yac_id, "2000-01-01T00:00:00", "2000-01-01T00:00:10")
    CALL yac_fdef_comp(yac_id, "comp", comp_id)
    CALL yac_fdef_grid("grid1", n_corners, cyclic, &
         & corner_lon, corner_lat, grid_id_out)
    CALL yac_fdef_grid("grid2", n_corners, cyclic, &
         & corner_lon, corner_lat, grid_id_in)
    CALL yac_fdef_points(grid_id_out, n_cells, YAC_LOCATION_CELL, &
         & cell_lon, cell_lat, point_id_out)
    CALL yac_fdef_points(grid_id_in, n_cells, YAC_LOCATION_CELL, &
         & cell_lon, cell_lat, point_id_in)

    CALL yac_fdef_field( &
      "out", comp_id, [point_id_out], 1, 1, "1", &
      YAC_TIME_UNIT_SECOND, field_id_out)
    CALL yac_fdef_field( &
      "in", comp_id, [point_id_in], 1, 1, "1", &
      YAC_TIME_UNIT_SECOND, field_id_in)

    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_average( &
      interp_stack, YAC_AVG_ARITHMETIC, 1)
    CALL yac_fdef_couple(yac_id, &
         "comp", "grid1", "out", &
         "comp", "grid2", "in", &
         "1", YAC_TIME_UNIT_SECOND, &
         YAC_REDUCTION_TIME_NONE, interp_stack, &
         src_lag = 0, tgt_lag = 0, &
         weight_file = weight_file_name, &
         mapping_side = 1, &
         scale_factor = 1.0d0, scale_summand = 0.0d0)
    CALL yac_ffree_interp_stack_config(interp_stack)

    CALL yac_fenddef(yac_id)
    CALL yac_ffinalize(yac_id)
  END SUBROUTINE gen_weight_file

  SUBROUTINE comp1_main()

    IMPLICIT NONE

    INTEGER :: field_out, avg_interp

    CALL yac_finit()
    ! this component defines the datetime
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fdef_datetime("2000-01-01T00:00:00", "2000-01-01T00:00:10")

    CALL yac_fset_grid_output_file("grid1", grid_output_filename)

    ! test yac_fpredef_comp:
    CALL yac_fpredef_comp("comp1", comp_id)
    CALL yac_fdef_comp_dummy()

    CALL yac_fdef_grid("grid1", n_corners, cyclic, &
         & corner_lon, corner_lat, grid_id)
    CALL yac_fdef_points(grid_id, n_cells, YAC_LOCATION_CELL, &
         & cell_lon, cell_lat, points_id)

    CALL yac_fdef_field( &
      "A", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_1)
    CALL yac_fdef_field( &
      "B", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_2)
    CALL yac_fdef_field( &
      "C", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_3)
    CALL yac_fdef_field( &
      "D", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_4)

    CALL yac_fdef_field( &
      "OUT", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_out)

    CALL yac_fget_interp_stack_config(avg_interp)
    CALL yac_fadd_interp_stack_config_average( &
      avg_interp, YAC_AVG_ARITHMETIC, 1)
    CALL yac_fdef_couple( &
         "comp2", "grid2", "4", &
         "comp1", "grid1", "D", &
         "1", YAC_TIME_UNIT_SECOND, &
         YAC_REDUCTION_TIME_ACCUMULATE, avg_interp, &
         src_lag = 0, tgt_lag = 0, &
         weight_file = weight_file_name_dummy, mapping_side = 1)
    CALL yac_ffree_interp_stack_config(avg_interp)

    CALL yac_fset_config_output_file( &
      test_config_filename_yaml, YAC_CONFIG_OUTPUT_FORMAT_YAML, &
      YAC_CONFIG_OUTPUT_SYNC_LOC_SYNC_DEF)
    CALL yac_fset_config_output_file( &
      test_config_filename_json, YAC_CONFIG_OUTPUT_FORMAT_JSON, &
      YAC_CONFIG_OUTPUT_SYNC_LOC_ENDDEF)

    CALL yac_fenddef( )

    DO WHILE (.TRUE.)
       timestep_counter = timestep_counter + 1
       data = 42.
       CALL yac_fput(field_1, SIZE(data), 1, 1, RESHAPE(data, [SIZE(data), 1, 1]), info, ierror)
       WRITE (0, *) "time of field_1: ", yac_fget_field_datetime(field_1)
       data = 43.
       CALL yac_fput(field_2, SIZE(data), 1, 1, RESHAPE(data, [SIZE(data), 1, 1]), info, ierror)
       data = 44.
       CALL yac_fput(field_3, SIZE(data), 1, 1, RESHAPE(data, [SIZE(data), 1, 1]), info, ierror)
       data = -1.
       CALL yac_fget(field_4, SIZE(data), 1, data, info, ierror)
       IF ((data(1,1) - 40.) > 1E-14) STOP 1
       data = 45.
       CALL yac_fput(field_out, SIZE(data), 1, 1, RESHAPE(data, [SIZE(data), 1, 1]), info, ierror)
       IF ( info .EQ. YAC_ACTION_PUT_FOR_RESTART) EXIT
    END DO

    CALL yac_ffinalize()

  END SUBROUTINE comp1_main

  SUBROUTINE comp2_main()
    IMPLICIT NONE
    INTEGER :: yac_id, avg_interp, i, field_in(9)
    CHARACTER (LEN=:), ALLOCATABLE :: config

    CALL yac_finit(yac_id)
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fpredef_comp(yac_id, "comp2", comp_id)
    CALL yac_fdef_comp_dummy(yac_id)

    CALL yac_fsync_def(yac_id)

    ! define fields after sync_def to tests the sync mechanism
    CALL yac_fdef_grid("grid2", n_corners, cyclic, &
         & corner_lon, corner_lat, grid_id)
    CALL yac_fdef_points(grid_id, n_cells, YAC_LOCATION_CELL, &
         & cell_lon, cell_lat, points_id)

    CALL yac_fset_grid_output_file(yac_id, "grid2", grid_output_filename)

    CALL yac_fdef_field( &
      "1", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_1)
    CALL yac_fdef_field( &
      "2", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_2)
    CALL yac_fdef_field( &
      "3", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_3)
    CALL yac_fdef_field( &
      "4", comp_id, [points_id], 1, 1, "1", YAC_TIME_UNIT_SECOND, field_4)

    CALL yac_fget_interp_stack_config(avg_interp)
    CALL yac_fadd_interp_stack_config_average( &
      avg_interp, YAC_AVG_ARITHMETIC, 1)

    CALL yac_fdef_couple(yac_id, &
         "comp1", "grid1", "A", &
         "comp2", "grid2", "1", &
         "1", YAC_TIME_UNIT_SECOND, &
         YAC_REDUCTION_TIME_ACCUMULATE, avg_interp)
    CALL yac_fdef_couple(yac_id, &
         "comp1", "grid1", "B", &
         "comp2", "grid2", "2", &
         "1", YAC_TIME_UNIT_SECOND, &
         YAC_REDUCTION_TIME_ACCUMULATE, avg_interp)
    CALL yac_fdef_couple(yac_id, &
         "comp1", "grid1", "C", &
         "comp2", "grid2", "3", &
         "1", YAC_TIME_UNIT_SECOND, &
         YAC_REDUCTION_TIME_ACCUMULATE, avg_interp)

    CALL yac_ffree_interp_stack_config(avg_interp)

    ! check all interpolation methods
    field_in(1) = gen_field_avg(yac_id)
    field_in(2) = gen_field_ncc(yac_id)
    field_in(3) = gen_field_nnn(yac_id)
    field_in(4) = gen_field_conserv(yac_id)
    field_in(5) = gen_field_spmap(yac_id)
    field_in(6) = gen_field_hcsbb(yac_id)
    field_in(7) = gen_field_user_file(yac_id)
    field_in(8) = gen_field_fixed(yac_id)
    field_in(9) = gen_field_creep(yac_id)

    CALL yac_fenddef( &
      yac_id, YAC_YAML_EMITTER_DEFAULT_F, config)
    DEALLOCATE(config)

    DO WHILE (.TRUE.)
       timestep_counter = timestep_counter + 1
       data = -1.0
       CALL yac_fget(field_1, SIZE(data), 1, data, info, ierror)
       WRITE (0, *) "time of field_1: ", yac_fget_field_datetime(field_1)
       IF ((data(1,1) - 42.) > 1E-14) STOP 1
       data = -1.0
       CALL yac_fget(field_2, SIZE(data), 1, data, info, ierror)
       IF ((data(1,1) - 43.) > 1E-14) STOP 1
       data = -1.0
       CALL yac_fget(field_3, SIZE(data), 1, data, info, ierror)
       IF ((data(1,1) - 44.) > 1E-14) STOP 1
       data = 40.
       CALL yac_fput(field_4, SIZE(data), 1, 1, RESHAPE(data, [SIZE(data), 1, 1]), info, ierror)
       DO i = 1, SIZE(field_in)
        CALL yac_fget(field_in(i), SIZE(data), 1, data, info, ierror)
       END DO
       IF ( info .EQ. YAC_ACTION_GET_FOR_RESTART) EXIT
    END DO

    CALL yac_ffinalize(yac_id)

  END SUBROUTINE comp2_main

  FUNCTION gen_field(yac_id, field_name, interp_stack) RESULT(field_id)
    INTEGER, INTENT(IN) :: yac_id
    CHARACTER(len=*), INTENT(IN) :: field_name
    INTEGER, INTENT(IN) :: interp_stack
    INTEGER :: field_id

    CALL yac_fdef_field( &
      field_name, comp_id, [points_id], 1, 1, "1", &
      YAC_TIME_UNIT_SECOND, field_id)
    CALL yac_fdef_couple( &
      yac_id, "comp1", "grid1", "OUT", "comp2", "grid2", field_name, &
      "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, interp_stack)
    CALL yac_ffree_interp_stack_config(interp_stack)
  END FUNCTION gen_field

  FUNCTION gen_field_avg(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_avg
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_average( &
      interp_stack, YAC_AVG_ARITHMETIC, 1)
    gen_field_avg = gen_field(yac_id, "avg", interp_stack)
  END FUNCTION gen_field_avg

  FUNCTION gen_field_ncc(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_ncc
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_ncc( &
      interp_stack, YAC_NCC_AVG, 1)
    gen_field_ncc = gen_field(yac_id, "ncc", interp_stack)
  END FUNCTION gen_field_ncc

  FUNCTION gen_field_nnn(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_nnn
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_nnn( &
      interp_stack, YAC_NNN_AVG, 3, 0.0d0, 1.0d0)
    gen_field_nnn = gen_field(yac_id, "nnn", interp_stack)
  END FUNCTION gen_field_nnn

  FUNCTION gen_field_conserv(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_conserv
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_conservative( &
      interp_stack, 1, 1, 1, YAC_CONSERV_DESTAREA)
    gen_field_conserv = gen_field(yac_id, "conserv", interp_stack)
  END FUNCTION gen_field_conserv

  FUNCTION gen_field_spmap(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_spmap
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_spmap( &
      interp_stack, 0.0d0, 0.0d0, YAC_SPMAP_AVG, YAC_SPMAP_NONE, 1.0d0, 1.0d0)
    gen_field_spmap = gen_field(yac_id, "spmap", interp_stack)
  END FUNCTION gen_field_spmap

  FUNCTION gen_field_hcsbb(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_hcsbb
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_hcsbb(interp_stack)
    gen_field_hcsbb = gen_field(yac_id, "hcsbb", interp_stack)
  END FUNCTION gen_field_hcsbb

  FUNCTION gen_field_user_file(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_user_file
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_user_file( &
      interp_stack, weight_file_name)
    gen_field_user_file = gen_field(yac_id, "user_file", interp_stack)
  END FUNCTION gen_field_user_file

  FUNCTION gen_field_fixed(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_fixed
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_fixed( &
      interp_stack, -1.0d0)
    gen_field_fixed = gen_field(yac_id, "fixed", interp_stack)
  END FUNCTION gen_field_fixed

  FUNCTION gen_field_creep(yac_id)
    INTEGER, INTENT(IN) :: yac_id
    INTEGER :: gen_field_creep
    INTEGER :: interp_stack
    CALL yac_fget_interp_stack_config(interp_stack)
    CALL yac_fadd_interp_stack_config_creep( &
      interp_stack, 0)
    gen_field_creep = gen_field(yac_id, "creep", interp_stack)
  END FUNCTION gen_field_creep

  SUBROUTINE error_exit ()

    USE mpi, ONLY : mpi_abort, MPI_COMM_WORLD
    USE utest

    INTEGER :: ierror

    CALL test ( .FALSE. )
    CALL stop_test
    CALL exit_tests
    CALL mpi_abort ( MPI_COMM_WORLD, 999, ierror )

  END SUBROUTINE error_exit

END PROGRAM main
