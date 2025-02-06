! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#include "test_macros.inc"

PROGRAM test_query_routines

  USE utest
  USE mpi
  use yaxt
  USE yac
  IMPLICIT NONE

  INTEGER :: global_rank, global_size, rank_type, ierror, i, j, rank, idx
  INTEGER :: yac_id, comp_id, comp_id_instance, grid_id, point_id, &
       field_id, field_id_instance
  INTEGER :: ref_nbr_comps, ref_nbr_grids
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: comp_name
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: ref_comp_name, ref_grid_names(3), &
                                     ref_field_names(2)
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: grid_name
  CHARACTER (LEN=YAC_MAX_CHARLEN) :: field_name
  TYPE(yac_string), ALLOCATABLE :: comp_names(:), comp_names_instance(:)
  TYPE(yac_string), ALLOCATABLE :: grid_names(:), grid_names_instance(:)
  TYPE(yac_string), ALLOCATABLE :: field_names(:), field_names_instance(:)
  TYPE(yac_string), ALLOCATABLE :: comp_grid_names(:)

  ! this avoids some compiler warnings
  allocate(comp_names(0), comp_names_instance(0), &
           grid_names(0), grid_names_instance(0), &
           field_names(0), field_names_instance(0), &
           comp_grid_names(0))

  CALL start_test("test_query_routines")

  CALL MPI_Init(ierror)
  CALL xt_initialize(MPI_COMM_WORLD)

  CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
  CALL yac_finit()
  CALL yac_finit(yac_id)

  CALL MPI_Comm_rank(MPI_COMM_WORLD, global_rank, ierror)
  CALL MPI_Comm_size(MPI_COMM_WORLD, global_size, ierror)


  !--------------------------------------------------
  ! definition of components, grid, points and fields
  !--------------------------------------------------

  rank_type = MOD(global_rank, 5)

  IF (rank_type == 0) THEN

    ! ranks of type 0 define no data
    CALL yac_fdef_comp_dummy()
    CALL yac_fdef_comp_dummy(yac_id)

  ELSE

    ! ranks of type 1 or higher define a component
    WRITE (comp_name , "('comp_',I0)") global_rank
    CALL yac_fdef_comp(comp_name, comp_id)
    CALL yac_fdef_comp(yac_id, comp_name, comp_id_instance)

    ! ranks of type 2 or higher define a grid and attach some
    ! meta data to their component
    IF (rank_type > 1) THEN

      CALL yac_fdef_component_metadata( &
        comp_name, TRIM(comp_name) // " METADATA_C")
      CALL yac_fdef_component_metadata( &
        yac_id, comp_name, TRIM(comp_name) // " METADATA_C_instance")

      WRITE (grid_name , "('grid_',I0,'_0')") global_rank
      CALL yac_fdef_grid( &
        grid_name, [2, 2], [0,0], [0.,180.], [-45.,45.], grid_id)

      CALL yac_fdef_points( &
        grid_id, [1,1], YAC_LOCATION_CELL, [90.], [0.], point_id)

      ! ranks of type 3 or higher define a second grid, attach some
      ! meta data to it
      IF (rank_type > 2) THEN

        WRITE (grid_name , "('grid_',I0,'_1')") global_rank
        CALL yac_fdef_grid( &
          grid_name, [2, 2], [0,0], [0.,180.], [-45.,45.], grid_id)

        CALL yac_fdef_points( &
          grid_id, [1,1], YAC_LOCATION_CELL, [90.], [0.], point_id)

        CALL yac_fdef_grid_metadata( &
          grid_name, TRIM(grid_name) // " METADATA_G")
        CALL yac_fdef_grid_metadata( &
          yac_id, grid_name, TRIM(grid_name) // " METADATA_G_instance")


        ! ranks of type 4 or higher define a third grid, attach some
        ! meta data to it and define two fields (one with meta data
        ! and one without)
        IF (rank_type > 3) THEN

          WRITE (grid_name , "('grid_',I0,'_2')") global_rank
          CALL yac_fdef_grid( &
            grid_name, [2, 2], [0,0], [0.,180.], [-45.,45.], grid_id)

          CALL yac_fdef_points( &
            grid_id, [1,1], YAC_LOCATION_CELL, [90.], [0.], point_id)

          CALL yac_fdef_grid_metadata( &
            grid_name, TRIM(grid_name) // " METADATA_G")
          CALL yac_fdef_grid_metadata( &
            yac_id, grid_name, TRIM(grid_name) // " METADATA_G_instance")

          WRITE (field_name , "('field_',I0,'_0')") global_rank
          CALL yac_fdef_field( &
            field_name, comp_id, [point_id], 1, 1, "PT5M", &
            YAC_TIME_UNIT_ISO_FORMAT, field_id)
          CALL yac_fdef_field( &
            field_name, comp_id_instance, [point_id], 1, 1, "PT5M", &
            YAC_TIME_UNIT_ISO_FORMAT, field_id_instance)

          WRITE (field_name , "('field_',I0,'_1')") global_rank
          CALL yac_fdef_field( &
            field_name, comp_id, [point_id], 1, 1, "PT5M", &
            YAC_TIME_UNIT_ISO_FORMAT, field_id)
          CALL yac_fdef_field( &
            field_name, comp_id_instance, [point_id], 1, 1, "PT5M", &
            YAC_TIME_UNIT_ISO_FORMAT, field_id_instance)
          CALL yac_fdef_field_metadata( &
            comp_name, grid_name, field_name, &
            TRIM(field_name) // " METADATA_F")
          CALL yac_fdef_field_metadata( &
            yac_id, comp_name, grid_name, field_name,&
            TRIM(field_name) // " METADATA_F_instance")

          !----------------------------------------
          ! test some field_id based query routines
          !----------------------------------------
          CALL test(field_id == yac_fget_field_id(comp_name, grid_name, field_name))
          CALL test(field_id_instance == yac_fget_field_id(yac_id, comp_name, grid_name, field_name))
          CALL test(comp_name == yac_fget_component_name(field_id))
          CALL test(comp_name == yac_fget_component_name(field_id_instance))
          CALL test(grid_name == yac_fget_grid_name(field_id))
          CALL test(grid_name == yac_fget_grid_name(field_id_instance))
          CALL test(field_name == yac_fget_field_name(field_id))
          CALL test(field_name == yac_fget_field_name(field_id_instance))
          CALL test("PT5M" == yac_fget_field_timestep(field_id))
          CALL test("PT5M" == yac_fget_field_timestep(field_id_instance))
          CALL test(1 == yac_fget_field_collection_size(field_id))
          CALL test(1 == yac_fget_field_collection_size(field_id_instance))
          CALL test(YAC_EXCHANGE_TYPE_NONE == yac_fget_field_role(field_id))
          CALL test(YAC_EXCHANGE_TYPE_NONE == yac_fget_field_role(field_id_instance))
        END IF
      END IF
    END IF
  END IF

  ! end definition phase and synchronise definitions across all processes
  CALL yac_fsync_def()
  CALL yac_fsync_def(yac_id)

  !--------------------
  ! test query routines
  !--------------------

  ref_nbr_comps = 4 * ((global_size - 1) / 5) + MOD(global_size - 1, 5)
  IF (ALLOCATED(comp_names)) DEALLOCATE(comp_names)
  comp_names = yac_fget_comp_names()
  IF (ALLOCATED(comp_names_instance)) DEALLOCATE(comp_names_instance)
  comp_names_instance = yac_fget_comp_names(yac_id)
  CALL test(SIZE(comp_names) == ref_nbr_comps)
  CALL test(SIZE(comp_names_instance) == ref_nbr_comps)

  ref_nbr_grids = &
    3 * ((global_size - 1) / 5) + MERGE(1, 0, MOD(global_size - 1, 5) > 2) + &
    MERGE(2, 0, MOD(global_size - 1, 5) > 3)

  IF (ALLOCATED(grid_names)) DEALLOCATE(grid_names)
  grid_names = yac_fget_grid_names()
  IF (ALLOCATED(grid_names_instance)) DEALLOCATE(grid_names_instance)
  grid_names_instance = yac_fget_grid_names(yac_id)
  CALL test(SIZE(grid_names) == ref_nbr_grids)
  CALL test(SIZE(grid_names_instance) == ref_nbr_grids)

  DO rank = 0, global_size - 1

    rank_type = MOD(rank, 5)
    WRITE (ref_comp_name, "('comp_',I0)") rank
    DO i = 1, 3
      WRITE (ref_grid_names(i), "('grid_',I0,'_',I0)") rank, i - 1
    END DO
    DO i = 1, 2
      WRITE (ref_field_names(i), "('field_',I0,'_',I0)") rank, i - 1
    END DO

    IF (rank_type == 0) THEN

      ! no component or grid is defined on this rank
      idx = 0
      DO i = 1, ref_nbr_comps
        IF (ref_comp_name == comp_names(i)%string) idx = i
      END DO
      CALL test(idx == 0)
      idx = 0
      DO i = 1, ref_nbr_comps
        IF (ref_comp_name == comp_names_instance(i)%string) idx = i
      END DO
      CALL test(idx == 0)
      DO i = 1, 3
        idx = 0
        DO j = 1, ref_nbr_grids
          IF (ref_grid_names(i) == grid_names(j)%string) idx = j
        END DO
        CALL test(idx == 0)
        idx = 0
        DO j = 1, ref_nbr_grids
          IF (ref_grid_names(i) == grid_names_instance(j)%string) idx = j
        END DO
        CALL test(idx == 0)
      END DO

    ELSE

      ! a component is defined on this process
      idx = 0
      DO i = 1, ref_nbr_comps
        IF (ref_comp_name == comp_names(i)%string) idx = i
      END DO
      CALL test(idx /= 0)
      idx = 0
      DO i = 1, ref_nbr_comps
        IF (ref_comp_name == comp_names_instance(i)%string) idx = i
      END DO
      CALL test(idx /= 0)

      IF (rank_type == 1) THEN

        ! the component of this rank has no meta data
        CALL test(.NOT. yac_fcomponent_has_metadata(ref_comp_name))
        CALL test(.NOT. yac_fcomponent_has_metadata(yac_id, ref_comp_name))

        ! no grid is registered in the coupling configuration on this rank
        DO i = 1, 3
          idx = 0
          DO j = 1, ref_nbr_grids
            IF (ref_grid_names(i) == grid_names(j)%string) idx = j
          END DO
          CALL test(idx == 0)
          idx = 0
          DO j = 1, ref_nbr_grids
            IF (ref_grid_names(i) == grid_names_instance(j)%string) idx = j
          END DO
          CALL test(idx == 0)
        END DO

        IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
        comp_grid_names = yac_fget_comp_grid_names(ref_comp_name)
        CALL test(SIZE(comp_grid_names) == 0)
        IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
        comp_grid_names = yac_fget_comp_grid_names(yac_id, ref_comp_name)
        CALL test(SIZE(comp_grid_names) == 0)

      ELSE

        ! the component has meta data on this rank
        CALL test(yac_fcomponent_has_metadata(ref_comp_name))
        CALL test(yac_fcomponent_has_metadata(yac_id, ref_comp_name))
        CALL test(yac_fget_component_metadata(ref_comp_name) == TRIM(ref_comp_name) // " METADATA_C")
        CALL test(yac_fget_component_metadata(yac_id, ref_comp_name) == TRIM(ref_comp_name) // " METADATA_C_instance")

        IF (rank_type == 2) THEN

          ! no grid is registered in the coupling configuration on this rank
          DO i = 1, 3
            idx = 0
            DO j = 1, ref_nbr_grids
              IF (ref_grid_names(i) == grid_names(j)%string) idx = j
            END DO
            CALL test(idx == 0)
            idx = 0
            DO j = 1, ref_nbr_grids
              IF (ref_grid_names(i) == grid_names_instance(j)%string) idx = j
            END DO
            CALL test(idx == 0)
          END DO

          IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
          comp_grid_names = yac_fget_comp_grid_names(ref_comp_name)
          CALL test(SIZE(comp_grid_names) == 0)
          IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
          comp_grid_names = yac_fget_comp_grid_names(yac_id, ref_comp_name)
          CALL test(SIZE(comp_grid_names) == 0)
        ELSE

          ! grid_*_0 is not registered in the coupling configuration
          idx = 0
          DO j = 1, ref_nbr_grids
            IF (ref_grid_names(1) == grid_names(j)%string) idx = j
          END DO
          CALL test(idx == 0)
          idx = 0
          DO j = 1, ref_nbr_grids
            IF (ref_grid_names(1) == grid_names_instance(j)%string) idx = j
          END DO
          CALL test(idx == 0)

          ! grid_*_1 is registered in the coupling configuration, has meta
          ! data, but no fields
          idx = 0
          DO j = 1, ref_nbr_grids
            IF (ref_grid_names(2) == grid_names(j)%string) idx = j
          END DO
          CALL test(idx /= 0)
          idx = 0
          DO j = 1, ref_nbr_grids
            IF (ref_grid_names(2) == grid_names_instance(j)%string) idx = j
          END DO
          CALL test(idx /= 0)
          CALL test(yac_fgrid_has_metadata(ref_grid_names(2)))
          CALL test(yac_fgrid_has_metadata(yac_id, ref_grid_names(2)))
          CALL test(yac_fget_grid_metadata(ref_grid_names(2)) == TRIM(ref_grid_names(2)) // " METADATA_G")
          CALL test(yac_fget_grid_metadata(yac_id, ref_grid_names(2)) == TRIM(ref_grid_names(2)) // " METADATA_G_instance")
          IF (ALLOCATED(field_names)) DEALLOCATE(field_names)
          field_names = yac_fget_field_names(ref_comp_name, ref_grid_names(2))
          CALL test(SIZE(field_names) == 0)
          IF (ALLOCATED(field_names_instance)) DEALLOCATE(field_names_instance)
          field_names_instance = &
            yac_fget_field_names(yac_id, ref_comp_name, ref_grid_names(2))
          CALL test(SIZE(field_names_instance) == 0)

          IF (rank_type == 3) THEN

            ! grid_*_2 is not defined on this rank
            idx = 0
            DO j = 1, ref_nbr_grids
              IF (ref_grid_names(3) == grid_names(j)%string) idx = j
            END DO
            CALL test(idx == 0)
            idx = 0
            DO j = 1, ref_nbr_grids
              IF (ref_grid_names(3) == grid_names_instance(j)%string) idx = j
            END DO
            CALL test(idx == 0)

            IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
            comp_grid_names = yac_fget_comp_grid_names(ref_comp_name)
            CALL test(SIZE(comp_grid_names) == 0)
            IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
            comp_grid_names = yac_fget_comp_grid_names(yac_id, ref_comp_name)
            CALL test(SIZE(comp_grid_names) == 0)

          ELSE

            ! grid_*_2 is registered in the coupling configuration, has
            ! meta data, and two fields
            idx = 0
            DO j = 1, ref_nbr_grids
              IF (ref_grid_names(3) == grid_names(j)%string) idx = j
            END DO
            CALL test(idx /= 0)
            idx = 0
            DO j = 1, ref_nbr_grids
              IF (ref_grid_names(3) == grid_names_instance(j)%string) idx = j
            END DO
            CALL test(idx /= 0)
            CALL test(yac_fgrid_has_metadata(ref_grid_names(3)))
            CALL test(yac_fgrid_has_metadata(yac_id, ref_grid_names(3)))
            CALL test(yac_fget_grid_metadata(ref_grid_names(3)) == TRIM(ref_grid_names(3)) // " METADATA_G")
            CALL test(yac_fget_grid_metadata(yac_id, ref_grid_names(3)) == TRIM(ref_grid_names(3)) // " METADATA_G_instance")
            IF (ALLOCATED(field_names)) DEALLOCATE(field_names)
            field_names = &
              yac_fget_field_names(ref_comp_name, ref_grid_names(3))
            CALL test(SIZE(field_names) == 2)
            IF (ALLOCATED(field_names_instance)) DEALLOCATE(field_names_instance)
            field_names_instance = &
              yac_fget_field_names(yac_id, ref_comp_name, ref_grid_names(3))
            CALL test(SIZE(field_names_instance) == 2)
            idx = 0
            DO i = 1, 2
              DO j = 1, 2
                IF (ref_field_names(i) == field_names(j)%string) idx = j
              END DO
              CALL test(idx /= 0)
              idx = 0
              DO j = 1, 2
                IF (ref_field_names(i) == field_names_instance(j)%string) idx = j
              END DO
              CALL test(idx /= 0)
            END DO
            CALL test(.NOT. yac_ffield_has_metadata(ref_comp_name, ref_grid_names(3), ref_field_names(1)))
            CALL test(.NOT. yac_ffield_has_metadata(yac_id, ref_comp_name, ref_grid_names(3), ref_field_names(1)))
            CALL test(yac_ffield_has_metadata(ref_comp_name, ref_grid_names(3), ref_field_names(2)))
            CALL test(yac_ffield_has_metadata(yac_id, ref_comp_name, ref_grid_names(3), ref_field_names(2)))
            CALL test(yac_fget_field_metadata(ref_comp_name, ref_grid_names(3), ref_field_names(2)) == TRIM(ref_field_names(2)) // " METADATA_F")
            CALL test(yac_fget_field_metadata(yac_id, ref_comp_name, ref_grid_names(3), ref_field_names(2)) == TRIM(ref_field_names(2)) // " METADATA_F_instance")

            IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
            comp_grid_names = yac_fget_comp_grid_names(ref_comp_name)
            CALL test(SIZE(comp_grid_names) == 1)
            CALL test(comp_grid_names(1)%string == ref_grid_names(3))
            IF (ALLOCATED(comp_grid_names)) DEALLOCATE(comp_grid_names)
            comp_grid_names = yac_fget_comp_grid_names(yac_id, ref_comp_name)
            CALL test(SIZE(comp_grid_names) == 1)
            CALL test(comp_grid_names(1)%string == ref_grid_names(3))

            CALL test(yac_fget_field_role(ref_comp_name, ref_grid_names(3), ref_field_names(1)) == YAC_EXCHANGE_TYPE_NONE)
            CALL test(yac_fget_field_role(yac_id, ref_comp_name, ref_grid_names(3), ref_field_names(1)) == YAC_EXCHANGE_TYPE_NONE)
            CALL test(yac_fget_field_collection_size(ref_comp_name, ref_grid_names(3), ref_field_names(1)) == 1)
            CALL test(yac_fget_field_collection_size(yac_id, ref_comp_name, ref_grid_names(3), ref_field_names(1)) == 1)
            CALL test(yac_fget_field_timestep(ref_comp_name, ref_grid_names(3), ref_field_names(1)) == "PT5M")
            CALL test(yac_fget_field_timestep(yac_id, ref_comp_name, ref_grid_names(3), ref_field_names(1)) == "PT5M")
          END IF
        END IF
      END IF
    END IF
  END DO

  CALL yac_ffinalize(yac_id)
  CALL yac_ffinalize()
  CALL xt_finalize()
  CALL MPI_Finalize(ierror)

  CALL stop_test

  CALL exit_tests

END PROGRAM test_query_routines
