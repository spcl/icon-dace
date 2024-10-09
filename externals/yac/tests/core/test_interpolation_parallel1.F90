! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

PROGRAM main

  USE utest
  USE yac_core
  USE mpi
  USE, INTRINSIC :: iso_c_binding

  IMPLICIT NONE

  INTERFACE

    ! fortran interface to C routine from dist_grid_utils
    FUNCTION yac_generate_basic_grid_reg2d_c( &
      name, coordinates_x, coordinates_y, num_cells, &
      local_start, local_count, with_halo) &
      BIND ( c, name='yac_generate_basic_grid_reg2d' )

      USE, INTRINSIC :: iso_c_binding

      CHARACTER(KIND=c_char)     :: name(*)
      REAL(kind=c_double)        :: coordinates_x(*)
      REAL(kind=c_double)        :: coordinates_y(*)
      INTEGER(kind=c_size_t)     :: num_cells(2)
      INTEGER(kind=c_size_t)     :: local_start(2)
      INTEGER(kind=c_size_t)     :: local_count(2)
      INTEGER(kind=c_int), value :: with_halo

      TYPE(c_ptr) :: yac_generate_basic_grid_reg2d_c
    END FUNCTION yac_generate_basic_grid_reg2d_c
  END INTERFACE

  CHARACTER(LEN=*), PARAMETER :: grid_name_src = "src_grid"
  CHARACTER(LEN=*), PARAMETER :: grid_name_tgt = "tgt_grid"
  INTEGER(kind=c_size_t), PARAMETER :: collection_size = 10

  INTEGER, PARAMETER :: max_opt_arg_len = 1024
  CHARACTER(max_opt_arg_len) :: reorder_type_arg
  INTEGER(KIND=c_int) :: reorder_type
  INTEGER :: arg_len

  INTEGER :: comm_rank, comm_size, ierror

  LOGICAL :: tgt_flag

  INTEGER :: split_comm

  REAL(KIND=c_double), PARAMETER :: &
    YAC_RAD = 0.01745329251994329576923690768489_c_double ! M_PI / 180

  CALL start_test('interpolation_parallel1')

  CALL yac_mpi_init_c()
  CALL yac_yaxt_init_c(MPI_COMM_WORLD)

  CALL test(COMMAND_ARGUMENT_COUNT() == 1)
  CALL GET_COMMAND_ARGUMENT(1, reorder_type_arg, arg_len)
  CALL test((reorder_type_arg == "src") .OR. (reorder_type_arg == "tgt"))
  reorder_type = &
    MERGE(YAC_MAPPING_ON_SRC, YAC_MAPPING_ON_TGT, reorder_type_arg == "src")

  CALL MPI_Comm_rank(MPI_COMM_WORLD, comm_rank, ierror)
  CALL MPI_Comm_size(MPI_COMM_WORLD, comm_size, ierror)
  CALL test(comm_size == 4)
  CALL MPI_Barrier(MPI_COMM_WORLD, ierror)

  ! split processes into source and target

  tgt_flag = comm_rank < 2

  CALL MPI_Comm_split( &
    MPI_COMM_WORLD, MERGE(1,0,tgt_flag), 0, split_comm, ierror)

  IF (tgt_flag) THEN
    CALL target_main(MPI_COMM_WORLD, split_comm, reorder_type)
  ELSE
    CALL source_main(MPI_COMM_WORLD, split_comm, reorder_type)
  END IF

  CALL MPI_Comm_free(split_comm, ierror)

  CALL yac_mpi_finalize_c()

  CALL stop_test
  CALL exit_tests

CONTAINS

  !
  ! The source grid is distributed among 2 processes
  !
  ! The global source grid has 5x4 cells:
  !
  !    24--44--25--45--26--46--27--47--28--48--29
  !    |       |       |       |       |       |
  !    34  15  36  16  38  17  40  18  42  19  43
  !    |       |       |       |       |       |
  !    18--33--19--35--20--37--21--39--22--41--23
  !    |       |       |       |       |       |
  !    23  10  25  11  27  12  29  13  31  14  32
  !    |       |       |       |       |       |
  !    12--22--13--24--14--26--15--28--16--30--17
  !    |       |       |       |       |       |
  !    12  05  14  06  16  07  18  08  20  09  21
  !    |       |       |       |       |       |
  !    06--11--07--13--08--15--09--17--10--19--11
  !    |       |       |       |       |       |
  !    01  00  03  01  05  02  07  03  09  04  10
  !    |       |       |       |       |       |
  !    00--01--01--02--02--04--03--06--04--08--05
  !
  SUBROUTINE source_main(global_comm, source_comm, reorder_type)

    INTEGER, INTENT(IN) :: global_comm
    INTEGER, INTENT(IN) :: source_comm
    INTEGER(KIND=c_int) :: reorder_type

    INTEGER :: my_source_rank, ierror
    INTEGER(KIND=c_size_t) :: i

    REAL(KIND=c_double), PARAMETER :: &
      coordinates_x(6) = (/0.0_c_double, 1.0_c_double, 2.0_c_double, &
                           3.0_c_double, 4.0_c_double, 5.0_c_double/) * YAC_RAD
    REAL(KIND=c_double), PARAMETER :: &
      coordinates_y(5) = (/0.0_c_double, 1.0_c_double, 2.0_c_double, &
                           3.0_c_double, 4.0_c_double/) * YAC_RAD
    INTEGER(KIND=c_size_t), PARAMETER :: num_cells(2) = (/5,4/)
    INTEGER(KIND=c_size_t), PARAMETER :: &
      local_start(2,2) = RESHAPE((/0,0,3,0/),(/2,2/))
    INTEGER(KIND=c_size_t), PARAMETER :: &
      local_count(2,2) = RESHAPE((/3,4,2,4/),(/2,2/))
    INTEGER(KIND=c_int), PARAMETER :: with_halo = 1

    INTEGER, PARAMETER :: &
      src_data_int(25,2) = &
        RESHAPE( &
          (/0,1,2,3,-1, &
            6,7,8,9,-1, &
            12,13,14,15,-1, &
            18,19,20,21,-1, &
            24,25,26,27,-1, &
            -1,3,4,5, &
            -1,9,10,11, &
            -1,15,16,17, &
            -1,21,22,23, &
            -1,27,28,29, &
            -1,-1,-1,-1,-1/), (/25,2/))
    REAL(KIND=c_double), ALLOCATABLE, TARGET :: src_data(:,:)
    INTEGER(KIND=c_size_t) :: data_size
    TYPE(c_ptr), TARGET :: src_field_(collection_size)
    TYPE(c_ptr), TARGET :: src_field_collection(collection_size)
    TYPE(c_ptr) :: src_field

    TYPE(c_ptr) :: src_grid, tgt_grid
    TYPE(c_ptr) :: grid_pair
    TYPE(c_ptr) :: interp_grid
    TYPE(c_ptr) :: interp_stack_config
    TYPE(c_ptr) :: interp_method_stack
    TYPE(c_ptr) :: interp_weights
    TYPE(c_ptr) :: interpolation

    CALL MPI_Comm_rank(source_comm, my_source_rank, ierror)

    src_grid = &
      yac_generate_basic_grid_reg2d_c( &
        TRIM(grid_name_src) // c_null_char, &
        coordinates_x, coordinates_y, num_cells, &
        local_start(:, my_source_rank + 1), &
        local_count(:, my_source_rank + 1), with_halo)
    tgt_grid = &
      yac_basic_grid_empty_new_c(TRIM(grid_name_tgt) // c_null_char)

    grid_pair = yac_dist_grid_pair_new_c(src_grid, tgt_grid, global_comm)

    interp_grid = &
      yac_interp_grid_new_c( &
           grid_pair, TRIM(grid_name_src) // c_null_char, &
           TRIM(grid_name_tgt) // c_null_char, 1_c_size_t, &
           (/YAC_LOC_CORNER/), (/-1_c_size_t/), (/-1_c_size_t/), &
           YAC_LOC_CORNER, -1_c_size_t, -1_c_size_t)

    interp_stack_config = yac_interp_stack_config_new_c()
    CALL yac_interp_stack_config_add_average_c( &
      interp_stack_config, YAC_INTERP_AVG_ARITHMETIC, 0_c_int)
    CALL yac_interp_stack_config_add_fixed_c( &
      interp_stack_config, 1337.0_c_double)

    interp_method_stack = &
      yac_interp_stack_config_generate_c(interp_stack_config)

    interp_weights = &
      yac_interp_method_do_search_c(interp_method_stack, interp_grid)

    interpolation = &
      yac_interp_weights_get_interpolation_c( &
        interp_weights, reorder_type, collection_size, &
        yac_interpolation_get_const_frac_mask_no_value_c(), &
        1.0_c_double, 0.0_c_double)

    ! -------------------
    ! set up source data
    ! -------------------

    ! src_data dimensions [collection_idx]
    !                     [pointset_idx]
    !                     [local_idx]
    data_size = yac_basic_grid_get_data_size_c(src_grid, YAC_LOC_CORNER)
    ALLOCATE(src_data(data_size, collection_size))
    DO i = 1, collection_size
      src_data(1:data_size,i) = &
        REAL( &
          src_data_int(1:data_size,my_source_rank+1) + INT(i * 30), c_double)
      src_field_(i) = c_loc(src_data(:,i))
      src_field_collection(i) = c_loc(src_field_(i))
    END DO
    src_field = c_loc(src_field_collection)

    CALL yac_interpolation_execute_put_c(interpolation, src_field)

    DEALLOCATE(src_data)
    CALL yac_interpolation_delete_c(interpolation)
    CALL yac_interp_weights_delete_c(interp_weights)
    CALL yac_interp_method_delete_c(interp_method_stack);
    CALL yac_interp_stack_config_delete_c(interp_stack_config)
    CALL yac_interp_grid_delete_c(interp_grid)
    CALL yac_dist_grid_pair_delete_c(grid_pair)
    CALL yac_basic_grid_delete_c(tgt_grid)
    CALL yac_basic_grid_delete_c(src_grid)

  END SUBROUTINE

  !
  ! The target grid is distributed among 2 processes
  !
  ! The global target grid has 6x3 cells:
  !
  !    21--39--22--40--23--41--24--42--25--43--26--44--27
  !    |       |       |       |       |       |       |
  !    27  12  29  13  31  14  33  15  35  16  37  17  38
  !    |       |       |       |       |       |       |
  !    14--26--15--28--16--30--17--32--18--34--19--36--20
  !    |       |       |       |       |       |       |
  !    14  06  16  07  18  08  20  09  22  10  24  11  25
  !    |       |       |       |       |       |       |
  !    07--13--08--15--09--17--10--19--11--21--12--23--13
  !    |       |       |       |       |       |       |
  !    01  00  03  01  05  02  07  03  09  04  11  05  12
  !    |       |       |       |       |       |       |
  !    00--00--01--02--02--04--03--06--04--08--05--10--06
  !
  SUBROUTINE target_main(global_comm, target_comm, reorder_type)


    INTEGER, INTENT(IN) :: global_comm
    INTEGER, INTENT(IN) :: target_comm
    INTEGER(KIND=c_int) :: reorder_type

    INTEGER :: my_target_rank, ierror
    INTEGER(KIND=c_size_t) :: i, j

    REAL(KIND=c_double), PARAMETER :: &
      coordinates_x(7) = (/0.5_c_double, 1.5_c_double, 2.5_c_double, &
                           3.5_c_double, 4.5_c_double, 5.5_c_double, &
                           6.5_c_double/) * YAC_RAD
    REAL(KIND=c_double), PARAMETER :: &
      coordinates_y(4) = (/0.5_c_double, 1.5_c_double, &
                           2.5_c_double, 3.5_c_double/) * YAC_RAD
    INTEGER(KIND=c_size_t), PARAMETER :: num_cells(2) = (/6,3/)
    INTEGER(KIND=c_size_t), PARAMETER :: &
      local_start(2,2) = RESHAPE((/0,0,3,0/),(/2,2/))
    INTEGER(KIND=c_size_t), PARAMETER :: &
      local_count(2,2) = RESHAPE((/3,3,3,3/),(/2,2/))
    INTEGER(KIND=c_int), PARAMETER :: with_halo = 0

    REAL(KIND=c_double), TARGET :: tgt_data(16, collection_size)
    REAL(KIND=c_double), PARAMETER :: &
      ref_tgt_data(16, 2) = &
        RESHAPE( &
          (/3.5_c_double,4.5_c_double,5.5_c_double,6.5_c_double, &
            9.5_c_double,10.5_c_double,11.5_c_double,12.5_c_double, &
            15.5_c_double,16.5_c_double,17.5_c_double,18.5_c_double, &
            21.5_c_double,22.5_c_double,23.5_c_double,24.5_c_double, &
            6.5_c_double,7.5_c_double,1337.0_c_double,1337.0_c_double, &
            12.5_c_double,13.5_c_double,1337.0_c_double,1337.0_c_double, &
            18.5_c_double,19.5_c_double,1337.0_c_double,1337.0_c_double, &
            24.5_c_double,25.5_c_double,1337.0_c_double,1337.0_c_double/), &
          (/16,2/))
    TYPE(c_ptr), TARGET :: tgt_field_collection(collection_size)
    TYPE(c_ptr) :: tgt_field
    REAL(kind=c_double) :: ref_tgt_value

    TYPE(c_ptr) :: src_grid, tgt_grid
    TYPE(c_ptr) :: grid_pair
    TYPE(c_ptr) :: interp_grid
    TYPE(c_ptr) :: interp_stack_config
    TYPE(c_ptr) :: interp_method_stack
    TYPE(c_ptr) :: interp_weights
    TYPE(c_ptr) :: interpolation

    CALL MPI_Comm_rank(target_comm, my_target_rank, ierror)

    tgt_grid = &
      yac_generate_basic_grid_reg2d_c( &
        TRIM(grid_name_tgt) // c_null_char, &
        coordinates_x, coordinates_y, num_cells, &
        local_start(:, my_target_rank + 1), &
        local_count(:, my_target_rank + 1), with_halo)
    src_grid = &
      yac_basic_grid_empty_new_c(TRIM(grid_name_src) // c_null_char)

    grid_pair = yac_dist_grid_pair_new_c(tgt_grid, src_grid, global_comm)

    interp_grid = &
      yac_interp_grid_new_c( &
           grid_pair, TRIM(grid_name_src) // c_null_char, &
           TRIM(grid_name_tgt) // c_null_char, 1_c_size_t, &
           (/YAC_LOC_CORNER/), (/-1_c_size_t/), (/-1_c_size_t/), &
           YAC_LOC_CORNER, -1_c_size_t, -1_c_size_t)

    interp_stack_config = yac_interp_stack_config_new_c()
    CALL yac_interp_stack_config_add_average_c( &
      interp_stack_config, YAC_INTERP_AVG_ARITHMETIC, 0_c_int)
    CALL yac_interp_stack_config_add_fixed_c( &
      interp_stack_config, 1337.0_c_double)

    interp_method_stack = &
      yac_interp_stack_config_generate_c(interp_stack_config)

    interp_weights = &
      yac_interp_method_do_search_c(interp_method_stack, interp_grid)

    interpolation = &
      yac_interp_weights_get_interpolation_c( &
        interp_weights, reorder_type, collection_size, &
        yac_interpolation_get_const_frac_mask_no_value_c(), &
        1.0_c_double, 0.0_c_double)

    tgt_data = -1.0_c_double
    DO i = 1, collection_size
      tgt_field_collection(i) = c_loc(tgt_data(:,i))
    END DO
    tgt_field = c_loc(tgt_field_collection)

    CALL yac_interpolation_execute_get_c(interpolation, tgt_field)

    DO i = 1, collection_size
      DO j = 1, 16_c_size_t
        IF (ref_tgt_data(j, my_target_rank+1) == 1337.0_c_double) THEN
          ref_tgt_value = 1337.0_c_double
        ELSE
          ref_tgt_value = &
            ref_tgt_data(j, my_target_rank+1) + &
            REAL(i * 30_c_size_t, c_double)
        END IF
        CALL test(tgt_data(j,i) == ref_tgt_value)
      END DO
    END DO

    CALL yac_interpolation_delete_c(interpolation)
    CALL yac_interp_weights_delete_c(interp_weights)
    CALL yac_interp_method_delete_c(interp_method_stack);
    CALL yac_interp_stack_config_delete_c(interp_stack_config)
    CALL yac_interp_grid_delete_c(interp_grid)
    CALL yac_dist_grid_pair_delete_c(grid_pair)
    CALL yac_basic_grid_delete_c(src_grid)
    CALL yac_basic_grid_delete_c(tgt_grid)

  END SUBROUTINE

END PROGRAM
