! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#define DUMMY_VALUE (-1337.0)

#include "test_macros.inc"

#ifndef TEST_PRECISION
#error "TEST_PRECISION is not defined"
#endif

#ifndef YAC_PTR_TYPE
#error "YAC_PTR_TYPE is not defined"
#endif

PROGRAM main

  USE utest
  USE yac
  USE yaxt
  USE mpi

  IMPLICIT NONE

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)   !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !< double precision
  INTEGER, PARAMETER :: wp = TEST_PRECISION              !< selected working precision

  ! M_PI / 180
  REAL(wp), PARAMETER :: YAC_RAD = 0.017453292519943295769_wp

  CHARACTER (LEN=YAC_MAX_CHARLEN), PARAMETER :: comp_A = 'comp_A'
  CHARACTER (LEN=YAC_MAX_CHARLEN), PARAMETER :: comp_B = 'comp_B'

  INTEGER :: rank, ierror

  INTEGER :: with_field_mask, exchange_type
  INTEGER :: group_comms(1)
  CHARACTER(len=YAC_MAX_CHARLEN) :: group_names(1)

  ! ===================================================================

  CALL MPI_Init(ierror)
  CALL xt_initialize(MPI_COMM_WORLD)
  CALL start_test('dummy_coupling6')

  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

  DO exchange_type = 1, 3
    DO with_field_mask = 1, 2

      SELECT CASE (rank)
        CASE (0)
          group_names(1) = 'yac'
          CALL yac_fmpi_handshake(MPI_COMM_WORLD, group_names, group_comms)
          CALL yac_finit_comm_dummy(group_comms(1))
          CALL MPI_Comm_free(group_comms(1), ierror)
        CASE (1)
          CALL run_comp_a(with_field_mask, exchange_type)
        CASE (2)
          CALL run_comp_b(with_field_mask, exchange_type)
        CASE DEFAULT
          CALL test(.FALSE.)
      END SELECT

    END DO
  END DO

  CALL xt_finalize()
  CALL MPI_Finalize(ierror)

  CALL stop_test
  CALL exit_tests

CONTAINS

  SUBROUTINE check_pair_comm(pair_comm, proc_name)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: pair_comm

    CHARACTER, INTENT(IN) :: proc_name

    CHARACTER(LEN=1) :: sendbuf
    CHARACTER(LEN=2) :: recvbuf
    INTEGER :: ierror

    sendbuf(1:1) = proc_name

    CALL MPI_Allgather( &
      sendbuf, 1, MPI_CHARACTER, recvbuf, 1, MPI_CHARACTER, pair_comm, ierror)

    CALL test(recvbuf == "AB")

  END SUBROUTINE

  SUBROUTINE run_comp_a(with_field_mask, exchange_type)

    USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_value, ieee_positive_inf

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: with_field_mask
    INTEGER, INTENT(IN) :: exchange_type

    INTEGER :: comp_id, grid_id, cell_point_id, vertex_point_id, &
               cell_mask_id, vertex_mask_id, field_id, multi_field_id, &
               dummy_field_id
    INTEGER :: pair_comm

    INTEGER, PARAMETER :: nbr_vertices = 4
    INTEGER, PARAMETER :: nbr_cells = 2
    INTEGER, DIMENSION(2), PARAMETER :: nbr_vertices_per_cell = (/3,3/)
    REAL(wp), DIMENSION(4), PARAMETER :: &
      x_vertices =                               &
        (/0.0 * YAC_RAD, -1.0 * YAC_RAD, 1.0 * YAC_RAD, 0.0 * YAC_RAD/)
    REAL(wp), DIMENSION(4), PARAMETER :: &
      y_vertices =                               &
        (/1.0 * YAC_RAD, 0.0 * YAC_RAD, 0.0 * YAC_RAD, -1.0 * YAC_RAD/)
    INTEGER, DIMENSION(6), PARAMETER :: cell_to_vertex = (/1,2,3, 2,4,3/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      x_cells = (/0.0 * YAC_RAD, 0.0 * YAC_RAD/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      y_cells = (/0.5 * YAC_RAD, -0.5 * YAC_RAD/)
    REAL(wp), TARGET :: cell_data(nbr_cells), vertex_data(nbr_vertices)
    REAL(wp), TARGET :: cell_frac_mask_data(nbr_cells), &
                        vertex_frac_mask_data(nbr_vertices)

    INTEGER, PARAMETER :: nbr_hor_points = 2
    DOUBLE PRECISION :: frac_mask_fallback_value = -1d0
    DOUBLE PRECISION :: frac_mask_fallback_value_multi
    INTEGER, PARAMETER :: collection_size = 1
    INTEGER, PARAMETER :: nbr_pointsets = 1
    INTEGER :: info, send_info, recv_info

    REAL(wp) :: dummy_cell_data(nbr_cells, collection_size)

    TYPE(YAC_PTR_TYPE) :: send_field_ptr(2, collection_size)
    TYPE(YAC_PTR_TYPE) :: frac_mask_ptr(2, collection_size)

    frac_mask_fallback_value_multi = &
      ieee_value(frac_mask_fallback_value_multi, ieee_positive_inf)

    ! initialise YAC default instance
    CALL yac_finit()
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fdef_datetime('2000-01-01T00:00:00', '2000-01-01T00:00:12')

    ! define local component
    CALL yac_fdef_comp(comp_A, comp_id)

    ! get communicator that contains both components
    CALL yac_fget_comps_comm((/comp_A, comp_B/), 2, pair_comm)

    ! check the pair_comm
    CALL check_pair_comm(pair_comm, 'A')
    CALL MPI_Comm_free(pair_comm, ierror)

    !   define local grid
    !
    !      0
    !     / \
    !    / o \
    !   /     \
    !  1-------2   Eq.
    !   \     /
    !    \ o /
    !     \ /
    !      3
    CALL yac_fdef_grid(                                              &
      'grid_A', nbr_vertices, nbr_cells, SUM(nbr_vertices_per_cell), &
      nbr_vertices_per_cell, x_vertices, y_vertices, cell_to_vertex, grid_id)

    ! center points in cells (needed e.g. for nearest neighbour)
    CALL yac_fdef_points( &
      grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, cell_point_id)
    ! vertex points
    CALL yac_fdef_points( &
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices, &
      vertex_point_id)

    ! masks
    CALL yac_fdef_mask( &
      grid_id, nbr_cells, YAC_LOCATION_CELL, (/.TRUE.,.FALSE./), cell_mask_id)
    CALL yac_fdef_mask( &
      grid_id, nbr_vertices, YAC_LOCATION_CORNER, (/1,1,1,1/), vertex_mask_id)

    ! define field
    IF (with_field_mask == 1) THEN
      CALL yac_fdef_field( &
        'A_to_B_src', comp_id, (/cell_point_id/), 1, &
        collection_size, "1", YAC_TIME_UNIT_SECOND, field_id)
      CALL yac_fdef_field( &
        'A_to_B_multi', comp_id, (/cell_point_id, vertex_point_id/), 2, &
        collection_size, "1", YAC_TIME_UNIT_SECOND, multi_field_id)
    ELSE
      CALL yac_fdef_field_mask( &
        'A_to_B_src', comp_id, (/cell_point_id/), (/cell_mask_id/), 1, &
        collection_size, "1", YAC_TIME_UNIT_SECOND, field_id)
      CALL yac_fdef_field_mask( &
        'A_to_B_multi', comp_id, (/cell_point_id, vertex_point_id/), &
        (/cell_mask_id, vertex_mask_id/), 2, collection_size, "1", &
        YAC_TIME_UNIT_SECOND, multi_field_id)
    END IF
    CALL yac_fenable_field_frac_mask( &
      comp_A, "grid_A", "A_to_B_src", frac_mask_fallback_value)
    CALL yac_fenable_field_frac_mask( &
      comp_A, "grid_A", "A_to_B_multi", frac_mask_fallback_value_multi)
    CALL yac_fdef_field( &
      'dummy_field', comp_id, (/cell_point_id/), 1, collection_size, "1", &
      YAC_TIME_UNIT_SECOND, dummy_field_id)

    ! generate coupling
    CALL yac_fenddef()

    CALL test(yac_fget_field_frac_mask_fallback_value(comp_A, "grid_A", "A_to_B_multi") == frac_mask_fallback_value_multi)

    ! move data from comp_A to comp_B
    cell_data = (/3.0, 4.0/)
    dummy_cell_data = RESHAPE((/-1.0, -1.0/), (/nbr_hor_points, collection_size/))
    vertex_data = (/3.0, 4.0, 5.0, 6.0/)
    cell_frac_mask_data = (/1.0, 1.0/)
    vertex_frac_mask_data = (/1.0, 1.0, 1.0, 1.0/)
    IF (exchange_type == 1) THEN
      CALL yac_fput(                                                            &
        field_id, nbr_hor_points, nbr_pointsets, collection_size,               &
        RESHAPE(cell_data, (/nbr_hor_points, nbr_pointsets, collection_size/)), &
        RESHAPE(                                                                &
          cell_frac_mask_data,                                                  &
          (/nbr_hor_points, nbr_pointsets, collection_size/)), info, ierror)
      CALL yac_fget(                                                      &
        dummy_field_id, nbr_hor_points, collection_size, dummy_cell_data, &
        info, ierror)
    ELSE IF (exchange_type == 2) THEN
      CALL yac_fput(                                                       &
        field_id, nbr_hor_points, collection_size,                         &
        RESHAPE(cell_data, (/nbr_hor_points, collection_size/)),           &
        RESHAPE(cell_frac_mask_data, (/nbr_hor_points, collection_size/)), &
        info, ierror)
      CALL yac_fget(                                                      &
        dummy_field_id, nbr_hor_points, collection_size, dummy_cell_data, &
        info, ierror)
    ELSE
      CALL yac_fexchange(                                                        &
        field_id, dummy_field_id, nbr_hor_points, nbr_pointsets, nbr_hor_points, &
        collection_size,                                                         &
        RESHAPE(cell_data, (/nbr_hor_points, nbr_pointsets, collection_size/)),  &
        RESHAPE(                                                                 &
          cell_frac_mask_data,                                                   &
          (/nbr_hor_points, nbr_pointsets, collection_size/)), dummy_cell_data,  &
        send_info, recv_info, ierror)
    END IF
    send_field_ptr(1,1)%p => cell_data
    send_field_ptr(2,1)%p => vertex_data
    frac_mask_ptr(1,1)%p => cell_frac_mask_data
    frac_mask_ptr(2,1)%p => vertex_frac_mask_data
    CALL yac_fput( &
      multi_field_id, 2, collection_size, send_field_ptr, frac_mask_ptr, &
      info, ierror)

    CALL yac_ffinalize()

  END SUBROUTINE run_comp_a

  SUBROUTINE run_comp_b(with_field_mask, exchange_type)

    USE, INTRINSIC :: ieee_arithmetic

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: with_field_mask
    INTEGER, INTENT(IN) :: exchange_type

    INTEGER :: default_instance_id, comp_id, grid_id_B, point_id_B, &
               field_id, multi_field_id
    INTEGER :: pair_comm

    INTEGER :: internal_comm
    INTEGER :: internal_instance_id, internal_comp_id, grid_id_C, point_id_C
    INTEGER :: field_id_src, field_id_tgt, dummy_field_id
    INTEGER :: interp_stack_nnn, interp_stack_fixed

    INTEGER, DIMENSION(2), PARAMETER :: nbr_vertices = (/2,3/)
    INTEGER, DIMENSION(2), PARAMETER :: cyclic = (/0,0/)
    INTEGER, DIMENSION(2), PARAMETER :: nbr_cells = (/1,2/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      x_vertices = (/-0.5 * YAC_RAD, 0.5 * YAC_RAD/)
    REAL(wp), DIMENSION(3), PARAMETER :: &
      y_vertices = (/-1.0 * YAC_RAD, 0.0 * YAC_RAD, 1.0 * YAC_RAD/)
    REAL(wp), DIMENSION(1), PARAMETER :: &
      x_cells = (/0.0 * YAC_RAD/)
    REAL(wp), DIMENSION(2), PARAMETER :: &
      y_cells = (/-0.5 * YAC_RAD, 0.5 * YAC_RAD/)

    REAL(wp), DIMENSION(6), PARAMETER :: &
      internal_x_vertices = &
        (/-0.4 * YAC_RAD, 0.4 * YAC_RAD, &
          -0.4 * YAC_RAD, 0.4 * YAC_RAD, &
          -0.4 * YAC_RAD, 0.4 * YAC_RAD/)
    REAL(wp), DIMENSION(6), PARAMETER :: &
      internal_y_vertices = &
        (/-0.7 * YAC_RAD, -0.7 * YAC_RAD, &
          0.0 * YAC_RAD,  0.0 * YAC_RAD, &
          0.8 * YAC_RAD,  0.8 * YAC_RAD/)

    INTEGER, PARAMETER :: nbr_hor_points = 2
    DOUBLE PRECISION :: frac_mask_fallback_value
    DOUBLE PRECISION, PARAMETER :: scale_factor = 2.0
    DOUBLE PRECISION, PARAMETER :: scale_summand = 0.5
    INTEGER, PARAMETER :: collection_size = 1
    REAL(wp), DIMENSION(nbr_hor_points,collection_size), TARGET :: recv_field_B
    REAL(wp), DIMENSION(nbr_hor_points,collection_size), TARGET :: frac_mask
    REAL(wp), DIMENSION(nbr_hor_points,collection_size), TARGET :: dummy_field
    REAL(wp), DIMENSION(nbr_hor_points,collection_size), TARGET :: dummy_frac_mask
    REAL(wp), DIMENSION(6,collection_size), TARGET :: recv_field_C
    INTEGER :: info

    INTEGER :: send_info, recv_info

    TYPE(YAC_PTR_TYPE) :: send_field_ptr(1, collection_size)
    TYPE(YAC_PTR_TYPE) :: frac_mask_ptr(1, collection_size)
    TYPE(YAC_PTR_TYPE) :: recv_field_ptr(collection_size)
    TYPE(YAC_PTR_TYPE) :: dummy_field_ptr(1, collection_size)
    TYPE(YAC_PTR_TYPE) :: dummy_frac_mask_ptr(1, collection_size)

    frac_mask_fallback_value = &
      ieee_value(frac_mask_fallback_value, ieee_positive_inf)

    ! initialise YAC default instance
    CALL yac_finit(default_instance_id)
    CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
    CALL yac_fdef_datetime( &
      default_instance_id, '2000-01-01T00:00:00', '2000-01-01T00:00:12')

    ! define local component
    CALL yac_fdef_comp(default_instance_id, comp_B, comp_id)

    ! get communicator that contains both components
    CALL yac_fget_comps_comm( &
      default_instance_id, (/comp_A, comp_B/), 2, pair_comm)

    ! check the pair_comm
    CALL check_pair_comm(pair_comm, 'B')
    CALL MPI_Comm_free(pair_comm, ierror)

    !   define local grid
    !
    !   4-------5
    !   |       |
    !   |   o   |
    !   |       |
    !   2-------3
    !   |       |
    !   |   o   |
    !   |       |
    !   0-------1
    CALL yac_fdef_grid( &
      'grid_B', nbr_vertices, cyclic, x_vertices, y_vertices, grid_id_B)

    ! center points in cells (needed e.g. for nearest neighbour)
    CALL yac_fdef_points( &
      grid_id_B, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, point_id_B)

    ! define field
    CALL yac_fdef_field( &
      'A_to_B_tgt', comp_id, (/point_id_B/), 1, 1, "1", YAC_TIME_UNIT_SECOND, &
      field_id)
    CALL yac_fdef_field( &
      'A_to_B_multi', comp_id, (/point_id_B/), 1, 1, "1", YAC_TIME_UNIT_SECOND, &
      multi_field_id)
    CALL yac_fdef_field( &
      'dummy_field', comp_id, (/point_id_B/), 1, 1, "1", YAC_TIME_UNIT_SECOND, &
      dummy_field_id)

    ! define interpolation stacks
    CALL yac_fget_interp_stack_config_from_string_yaml( &
           "- nnn:" // ACHAR(10) // &
           "    n: 1", interp_stack_nnn)
    CALL yac_fget_interp_stack_config_from_string_json( &
           '[{"fixed": {"user_value": -2}}]', interp_stack_fixed)

    ! define couplings
    CALL yac_fdef_couple( &
      default_instance_id, &
      comp_A, "grid_A", "A_to_B_src", &
      comp_B, "grid_B", "A_to_B_tgt", &
      "PT01.000S", YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, &
      interp_stack_nnn)
    CALL yac_fdef_couple( &
      default_instance_id, &
      comp_A, "grid_A", "A_to_B_multi", &
      comp_B, "grid_B", "A_to_B_multi", &
      "PT01.000S", YAC_TIME_UNIT_ISO_FORMAT, YAC_REDUCTION_TIME_NONE, &
      interp_stack_fixed)
    CALL yac_ffree_interp_stack_config(interp_stack_nnn)
    CALL yac_ffree_interp_stack_config(interp_stack_fixed)

    ! generate coupling
    CALL yac_fenddef(default_instance_id)

    CALL test(yac_fget_field_frac_mask_fallback_value(default_instance_id, comp_A, "grid_A", "A_to_B_multi") == frac_mask_fallback_value)

    ! move data from comp_A to comp_B
    recv_field_B = DUMMY_VALUE
    recv_field_ptr(1)%p => recv_field_B(:,1)
    dummy_field = DUMMY_VALUE
    dummy_frac_mask = 1.0_wp
    CALL yac_fget( &
      multi_field_id, collection_size, recv_field_ptr, info, ierror)
    CALL test(ALL(recv_field_B == -2.0))
    IF (exchange_type == 1) THEN
      CALL yac_fexchange( &
        dummy_field_id, field_id, &
        nbr_hor_points, 1, nbr_hor_points, collection_size, &
        RESHAPE(dummy_field, (/nbr_hor_points, 1, collection_size/)), &
        recv_field_B, send_info, recv_info, ierror)
    ELSE
      CALL yac_fget( &
        field_id, nbr_hor_points, collection_size, recv_field_B, info, ierror)
    END IF

    ! do some additional tests with the dummy field
    CALL yac_fput(                                                  &
      dummy_field_id, nbr_hor_points, 1, collection_size,           &
      RESHAPE(dummy_field, (/nbr_hor_points, 1, collection_size/)), &
      info, ierror)
    CALL yac_fput(                                     &
      dummy_field_id, nbr_hor_points, collection_size, &
      dummy_field, info, ierror)
    CALL yac_fput(                                                      &
      dummy_field_id, nbr_hor_points, 1, collection_size,               &
      RESHAPE(dummy_field, (/nbr_hor_points, 1, collection_size/)),     &
      RESHAPE(dummy_frac_mask, (/nbr_hor_points, 1, collection_size/)), &
      info, ierror)
    CALL yac_fput(                                     &
      dummy_field_id, nbr_hor_points, collection_size, &
      dummy_field, dummy_frac_mask, info, ierror)
    dummy_field_ptr(1,1)%p => dummy_field(:,1)
    dummy_frac_mask_ptr(1,1)%p => dummy_frac_mask(:,1)
    CALL yac_fput(                        &
      dummy_field_id, 1, collection_size, &
      dummy_field_ptr,  info, ierror)
    CALL yac_fput(                          &
      dummy_field_id, 1, collection_size,   &
      dummy_field_ptr, dummy_frac_mask_ptr, &
      info, ierror)

    !-----------------------------------------------------
    ! run internal instance
    !-----------------------------------------------------

      CALL yac_fget_comp_comm(comp_id, internal_comm)

      ! initialise internal YAC instance
      CALL yac_finit_comm(internal_comm, internal_instance_id)
      CALL yac_fdef_datetime( &
        internal_instance_id, '2000-01-01T00:00:00', '2000-01-01T00:00:12')

      ! define internal component
      CALL yac_fdef_comp(internal_instance_id, comp_B, internal_comp_id)

      !   define local grid
      !
      !   4-------5
      !   |       |
      !   |   o   |
      !   |       |
      !   2-------3
      !   |       |
      !   |   o   |
      !   |       |
      !   0-------1
      CALL yac_fdef_grid(               &
        'grid_C', nbr_vertices(1) * nbr_vertices(2), &
        internal_x_vertices, internal_y_vertices, grid_id_C)

      ! vertices points (needed e.g. for nearest neighbour)
      CALL yac_fdef_points(                           &
        grid_id_C, nbr_vertices(1) * nbr_vertices(2), YAC_LOCATION_CORNER, &
        internal_x_vertices, internal_y_vertices, point_id_C)

      ! define field
      CALL yac_fdef_field( &
        'B_to_C', internal_comp_id, (/point_id_B/), 1, &
        collection_size, "1", YAC_TIME_UNIT_SECOND, field_id_src)
      CALL yac_fenable_field_frac_mask( &
        internal_instance_id, comp_B, "grid_B", "B_to_C", &
        frac_mask_fallback_value)
      CALL yac_fdef_field( &
        'B_to_C', internal_comp_id, (/point_id_C/), 1, collection_size, &
        "1", YAC_TIME_UNIT_SECOND, field_id_tgt)

      ! define interpolation stacks
      CALL yac_fget_interp_stack_config(interp_stack_nnn)
      CALL yac_fadd_interp_stack_config_nnn( &
        interp_stack_nnn, YAC_NNN_AVG, 1, 0.0_dp, 0.0_dp)

      ! define couplings
      CALL yac_fdef_couple( &
        internal_instance_id, &
        comp_B, "grid_B", "B_to_C", &
        comp_B, "grid_C", "B_to_C", &
        "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_NONE, interp_stack_nnn, &
        scale_factor=scale_factor, scale_summand=scale_summand)
      CALL yac_ffree_interp_stack_config(interp_stack_nnn)

      ! generate internal coupling
      CALL yac_fenddef(internal_instance_id)

      frac_mask = &
        RESHAPE((/1.0_wp, 0.0_wp/), (/nbr_hor_points, collection_size/))

      recv_field_C = DUMMY_VALUE

      IF (exchange_type == 1) THEN

        ! move data from grid B to grid C
        recv_field_C = DUMMY_VALUE
        CALL yac_fexchange( &
          field_id_src, field_id_tgt, nbr_hor_points, 6, collection_size,&
          RESHAPE(recv_field_B, (/nbr_hor_points, collection_size/)), &
          frac_mask, recv_field_C, send_info, recv_info, ierror)

        IF (with_field_mask == 1) THEN
          CALL test(ALL(RESHAPE(recv_field_C, (/6/)) == scale_factor * (/4.0_wp, 4.0_wp, 4.0_wp, 4.0_wp, REAL(frac_mask_fallback_value, wp), REAL(frac_mask_fallback_value, wp)/) + scale_summand))
        ELSE
          CALL test(ALL(RESHAPE(recv_field_C, (/6/)) == scale_factor * (/3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, REAL(frac_mask_fallback_value, wp), REAL(frac_mask_fallback_value, wp)/) + scale_summand))
        END IF

      ELSE

        send_field_ptr(1,1)%p => recv_field_B(:,1)
        frac_mask_ptr(1,1)%p => frac_mask(:,1)
        recv_field_ptr(1)%p => recv_field_C(:,1)

        ! move data from grid B to grid C
        CALL yac_fexchange( &
          field_id_src, field_id_tgt, 1, collection_size, &
          send_field_ptr, frac_mask_ptr, recv_field_ptr, &
          send_info, recv_info, ierror)

        IF (with_field_mask == 1) THEN
          CALL test(ALL(RESHAPE(recv_field_C, (/6/)) == scale_factor * (/4.0_wp, 4.0_wp, 4.0_wp, 4.0_wp, REAL(frac_mask_fallback_value, wp), REAL(frac_mask_fallback_value, wp)/) + scale_summand))
        ELSE
          CALL test(ALL(RESHAPE(recv_field_C, (/6/)) == scale_factor * (/3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, REAL(frac_mask_fallback_value, wp), REAL(frac_mask_fallback_value, wp)/) + scale_summand))
        END IF

      END IF

      CALL yac_fcleanup(internal_instance_id)

      CALL MPI_Comm_free(internal_comm, ierror)

    !-----------------------------------------------------

    CALL yac_ffinalize(default_instance_id)

  END SUBROUTINE

END PROGRAM
