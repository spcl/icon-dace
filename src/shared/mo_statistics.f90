! Set of methods for simple statistics
! NOTE: in order to get correct results make sure you provide the proper in_subset!
!
!
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

!----------------------------
#include "omp_definitions.inc"
!----------------------------
MODULE mo_statistics
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp,sp
  USE mo_exception,          ONLY: warning, finish
#ifdef _OPENMP
  USE omp_lib
#endif

  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed
  USE mo_mpi,                ONLY: process_mpi_stdio_id, get_my_mpi_work_communicator, p_max, p_min, &
    & my_process_is_mpi_parallel, p_sum, my_process_is_mpi_seq
  !   USE mo_io_units,           ONLY: nnml, filename_max
  !   USE mo_namelist,           ONLY: position_nml, open_nml, positioned
  USE mo_impl_constants, ONLY: on_cells, on_edges, on_vertices
  USE mo_math_types,     ONLY: t_geographical_coordinates
  USE mo_math_constants, ONLY: rad2deg
#ifdef _OPENACC
  USE mo_mpi,            ONLY: i_am_accel_node
#endif
  USE mo_fortran_tools,  ONLY: set_acc_host_or_device

  IMPLICIT NONE

  PRIVATE
#define VerticalDim_Position 2

  !-------------------------------------------------------------------------
  ! NOTE: in order to get correct results make sure you provide the proper in_subset (ie, owned)!
  PUBLIC :: global_minmaxmean, subset_sum, add_fields_3d
  PUBLIC :: add_fields
  PUBLIC :: L2Norm, LInfNorm
  PUBLIC :: accumulate_mean, levels_horizontal_mean, horizontal_mean, total_mean
  PUBLIC :: horizontal_sum
  PUBLIC :: print_value_location
  PUBLIC :: add_verticallyIntegrated_field, verticallyIntegrated_field
  PUBLIC :: add_verticalSum_field
  PUBLIC :: gather_sums

  ! simple min max mean (no weights)
  ! uses the range in_subset
  ! used from the debug print
  INTERFACE global_minmaxmean
    MODULE PROCEDURE MinMaxMean_2D
    MODULE PROCEDURE MinMaxMean_3D_AllLevels
    MODULE PROCEDURE MinMaxMean_2D_InRange
    MODULE PROCEDURE MinMaxMean_3D_AllLevels_InRange
  END INTERFACE global_minmaxmean

  INTERFACE L2Norm
    MODULE PROCEDURE L2Norm_2D_InRange
  END INTERFACE L2Norm

  INTERFACE LInfNorm
    MODULE PROCEDURE LInfNorm_2D_InRange
  END INTERFACE LInfNorm

  INTERFACE verticallyIntegrated_field
    MODULE PROCEDURE verticallyIntegrated_field_weighted
    MODULE PROCEDURE verticallyIntegrated_field_notWeighted
  END INTERFACE verticallyIntegrated_field

  ! weighted total sum, uses indexed subset
  ! used for calculating total fluxes acros given paths
  INTERFACE subset_sum
    MODULE PROCEDURE Sum_3D_AllLevels_3Dweights_InIndexed
    MODULE PROCEDURE Sum_2D_2Dweights_InRange
    MODULE PROCEDURE Sum_2D_InRange
    ! MODULE PROCEDURE globalspace_3d_sum_max_level_array
  END INTERFACE subset_sum

  ! function, restuns the total weighted mean across all levels
  INTERFACE total_mean
    MODULE PROCEDURE TotalWeightedMean_3D_InRange_3Dweights
  END INTERFACE total_mean

  ! for each level, adds to the level accumulated_mean the weighted spatial level mean
  INTERFACE accumulate_mean
    MODULE PROCEDURE AccumulateMean_3D_EachLevel_InRange_2Dweights
  END INTERFACE accumulate_mean

  ! for each level, gives the weighted spatial level mean
  INTERFACE levels_horizontal_mean
    MODULE PROCEDURE LevelHorizontalMean_3D_InRange_2Dweights
    MODULE PROCEDURE LevelHorizontalMean_3D_InRange_3Dweights
    MODULE PROCEDURE HorizontalMean_2D_InRange_2Dweights
  END INTERFACE levels_horizontal_mean

  ! the same as above, but better name
  INTERFACE horizontal_mean
    MODULE PROCEDURE LevelHorizontalMean_3D_InRange_2Dweights
    MODULE PROCEDURE LevelHorizontalMean_3D_InRange_3Dweights
    MODULE PROCEDURE HorizontalMean_2D_InRange_2Dweights
  END INTERFACE horizontal_mean

  INTERFACE horizontal_sum
    MODULE PROCEDURE LevelHorizontalSum_3D_InRange_2Dweights
    MODULE PROCEDURE LevelHorizontalSum_3D_InRange_3Dweights
  END INTERFACE horizontal_sum


  INTERFACE gather_sums
    MODULE PROCEDURE gather_sums_0D
    MODULE PROCEDURE gather_sums_1D
  END INTERFACE gather_sums

  PUBLIC :: construct_statistic_objects, destruct_statistic_objects
  PUBLIC :: new_statistic, delete_statistic
  PUBLIC :: add_statistic_to, add_data_to
!   PUBLIC :: MIN, MAX, mean
  PUBLIC :: max_statistic_of, min_statistic_of
  PUBLIC :: mean_statistic_of

  PUBLIC :: add_max_ratio

  PUBLIC :: new, delete

  PUBLIC :: t_statistic

  PUBLIC :: time_avg

  !-------------------------------------------------------------------------
  TYPE :: t_statistic
    INTEGER :: id
  END TYPE t_statistic
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! Parameters
  INTEGER, PARAMETER :: add_value = 1
  INTEGER, PARAMETER :: add_max_ratio = 2

  !--------------------------------------------------------------
  ! TYPE definitions
  !> Basic statistics type
  TYPE t_data_statistics
    LOGICAL :: is_active
    INTEGER :: no_of_values
    REAL(wp) :: sum_of_values
    REAL(wp) :: max_of_values
    REAL(wp) :: min_of_values
    ! distribution stitistics
    ! at the moment just a
    INTEGER :: no_of_bars
    REAL(wp) :: min_bars_value
    REAL(wp) :: max_bars_value
    INTEGER, POINTER :: no_of_values_in_bar(:)

    INTEGER :: mode ! defines how we insert the values,
    ! ADD_VALUE
    ! ADD_MAX_RATIO, input twos values

  END TYPE t_data_statistics
  TYPE t_data_statistics_2d
    LOGICAL :: is_active
    INTEGER :: no_of_values
    REAL(wp), POINTER :: sum_of_values(:,:)
    REAL(wp), POINTER :: max_of_values(:,:)
    REAL(wp), POINTER :: min_of_values(:,:)
    INTEGER :: no_of_bars
    REAL(wp) :: min_bars_value
    REAL(wp) :: max_bars_value
    INTEGER, POINTER :: no_of_values_in_bar(:)
    INTEGER :: mode
  END TYPE t_data_statistics_2d

  TYPE t_data_statistics_3d
    LOGICAL :: is_active
    INTEGER :: no_of_values
    REAL(wp), POINTER :: sum_of_values(:,:,:)
    REAL(wp), POINTER :: max_of_values(:,:,:)
    REAL(wp), POINTER :: min_of_values(:,:,:)
    INTEGER :: no_of_bars
    REAL(wp) :: min_bars_value
    REAL(wp) :: max_bars_value
    INTEGER, POINTER :: no_of_values_in_bar(:)
    INTEGER :: mode
  END TYPE t_data_statistics_3d

  TYPE t_data_statistics_collection
    INTEGER, ALLOCATABLE :: stats_1d(:)
    INTEGER, ALLOCATABLE :: stats_2d(:)
    INTEGER, ALLOCATABLE :: stats_3d(:)
  END TYPE

  !--------------------------------------------------------------
  !> The maximum number of statistic objects.
  INTEGER, PARAMETER ::  max_no_of_statistic_objects = 50
  !> The array of the statistic objects.
  TYPE(t_data_statistics),    ALLOCATABLE, TARGET :: statistic_object(:)
  TYPE(t_data_statistics_2d), ALLOCATABLE, TARGET :: statistic_object_2d(:)
  TYPE(t_data_statistics_3d), ALLOCATABLE, TARGET :: statistic_object_3d(:)
  !> The number of allocated statistic objects.
  INTEGER :: no_of_allocated_statistics = 0        ! the size of the statistic objects array
  !> The number of actual active statistic objects.
  INTEGER :: active_statistics
  !> The maximum id of the active statistics
  INTEGER :: max_active_statistics
  !> True if the statistic object is active.

  INTERFACE new
    MODULE PROCEDURE new_statistic_operator
  END INTERFACE
  INTERFACE delete
    MODULE PROCEDURE delete_statistic_operator
  END INTERFACE
  INTERFACE add_data_to
    MODULE PROCEDURE add_data_one_value_real
    MODULE PROCEDURE add_data_one_value_int
  END INTERFACE
!   INTERFACE MIN
!     MODULE PROCEDURE min_statistic
!   END INTERFACE
!   INTERFACE MAX
!     MODULE PROCEDURE max_statistic
!   END INTERFACE
!   INTERFACE mean
!     MODULE PROCEDURE mean_statistic
!   END INTERFACE
  INTERFACE new_statistic
    MODULE PROCEDURE new_statistic_no_bars
    MODULE PROCEDURE new_statistic_with_bars
    MODULE PROCEDURE new_statistic_with_pointers
  END INTERFACE
  INTERFACE init_statistics
    MODULE PROCEDURE init_statistic_object
    MODULE PROCEDURE init_statistic_object_with_pointers
  END INTERFACE
  INTERFACE add_statistic_to
    MODULE PROCEDURE add_statistic_one_value
    MODULE PROCEDURE add_statistic_two_values
  END INTERFACE

  INTERFACE add_fields
    MODULE PROCEDURE add_fields_3d
    MODULE PROCEDURE add_fields_2d
    MODULE PROCEDURE add_fields_3d_sp
    MODULE PROCEDURE add_fields_2d_sp
  END INTERFACE add_fields

  INTERFACE print_value_location
    MODULE PROCEDURE print_2Dvalue_location
    MODULE PROCEDURE print_3Dvalue_location
  END INTERFACE print_value_location

 CHARACTER(LEN=*), PARAMETER :: module_name="mo_statistics"

CONTAINS

  !-----------------------------------------------------------------------
  FUNCTION MinMaxMean_2D(values) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:)
    REAL(wp) :: minmaxmean(3)
    CALL warning(module_name, "not available without subset input")
    minmaxmean(:) = 1234567890
  END FUNCTION MinMaxMean_2D
  !-----------------------------------------------------------------------
  FUNCTION MinMaxMean_3D_AllLevels(values, start_level, end_level) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:,:)
    INTEGER, OPTIONAL :: start_level, end_level
    REAL(wp) :: minmaxmean(3)

    CALL warning(module_name, "not available without subset input")
    minmaxmean(:) = 1234567890

  END FUNCTION MinMaxMean_3D_AllLevels
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the min max mean in a 3D array in a given range subset
  ! The results are over all levels
  SUBROUTINE print_2Dvalue_location(values, seek_value, in_subset)
    REAL(wp), INTENT(in) :: values(:,:)
    REAL(wp), INTENT(in) :: seek_value
    TYPE(t_subset_range), TARGET, INTENT(in) :: in_subset

    INTEGER :: block, start_index, end_index, idx
    TYPE(t_geographical_coordinates), POINTER ::  geocoordinates(:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':print_cell_value_location'

!     IF (in_subset%no_of_holes > 0) CALL warning(module_name, "there are holes in the subset")

    ! get lon, lat
    SELECT CASE (in_subset%entity_location)
    CASE(on_cells)
      geocoordinates => in_subset%patch%cells%center
    CASE(on_edges)
      geocoordinates => in_subset%patch%edges%center
    CASE(on_vertices)
      geocoordinates => in_subset%patch%verts%vertex
    CASE default
      CALL finish(method_name, "unknown subset%entity_location")
    END SELECT

    DO block = in_subset%start_block, in_subset%end_block
      CALL get_index_range(in_subset, block, start_index, end_index)
      DO idx = start_index, end_index
        IF (values(idx, block) == seek_value) THEN
          WRITE(0,*) "Value ", seek_value, &
            & " found at lon=",  geocoordinates(idx, block)%lon * rad2deg, &
            & " lat=",  geocoordinates(idx, block)%lat * rad2deg
        ENDIF
      ENDDO
    ENDDO

  END SUBROUTINE print_2Dvalue_location
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  ! Returns the min max mean in a 3D array in a given range subset
  ! The results are over all levels
  SUBROUTINE print_3Dvalue_location(values, seek_value, in_subset, start_level, end_level)
    REAL(wp), INTENT(in) :: values(:,:,:)
    REAL(wp), INTENT(in) :: seek_value
    TYPE(t_subset_range), TARGET, INTENT(in) :: in_subset
    INTEGER, OPTIONAL, INTENT(in) :: start_level, end_level

    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    TYPE(t_geographical_coordinates), POINTER ::  geocoordinates(:,:)
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':print_cell_value_location'

!     IF (in_subset%no_of_holes > 0) CALL warning(module_name, "there are holes in the subset")

    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = end_level
    ELSE
      end_vertical = SIZE(values, VerticalDim_Position)
    ENDIF
    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")

    ! get lon, lat
    SELECT CASE (in_subset%entity_location)
    CASE(on_cells)
      geocoordinates => in_subset%patch%cells%center
    CASE(on_edges)
      geocoordinates => in_subset%patch%edges%center
    CASE(on_vertices)
      geocoordinates => in_subset%patch%verts%vertex
    CASE default
      CALL finish(method_name, "unknown subset%entity_location")
    END SELECT

    ! init the min, max values
    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = start_vertical, MIN(end_vertical, in_subset%vertical_levels(idx,block))
            IF (values(idx, level, block) == seek_value) THEN
              WRITE(0,*) "Value ", seek_value, &
                & " found at lon=",  geocoordinates(idx, block)%lon * rad2deg, &
                & " lat=",  geocoordinates(idx, block)%lat * rad2deg, &
                & " level=", level
           ENDIF
          ENDDO
        ENDDO
      ENDDO

    ELSE ! no in_subset%vertical_levels

      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = start_vertical, end_vertical
            IF (values(idx, level, block) == seek_value) THEN
              WRITE(0,*) "Value ", seek_value, &
                & " found at lon=",  geocoordinates(idx, block)%lon, &
                & " lat=",  geocoordinates(idx, block)%lat, &
                & " level=", level
            ENDIF
          ENDDO
        ENDDO
      ENDDO

    ENDIF

  END SUBROUTINE print_3Dvalue_location
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the l2 norm for a 2D array in a given range subset
  FUNCTION LInfNorm_2D_InRange(values, in_subset) result(lInfnorm)
    REAL(wp), INTENT(in) :: values(:,:)
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp) :: lInfnorm

    REAL(wp) :: minmaxmean(3)

    minmaxmean = global_minmaxmean(values, in_subset)
    lInfnorm = MAX(ABS(minmaxmean(1)), ABS(minmaxmean(2)))

  END FUNCTION LInfNorm_2D_InRange
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the l2 norm for a 2D array in a given range subset
  FUNCTION L2Norm_2D_InRange(values, in_subset) result(l2norm)
    REAL(wp), INTENT(in) :: values(:,:)
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp) :: l2norm

    REAL(wp) :: sumOfSquares
    INTEGER :: block, start_index, end_index, idx

    IF (in_subset%no_of_holes > 0) CALL warning(module_name, "there are holes in the subset")
    sumOfSquares = 0._wp

    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index, idx), reduction(+:sumOfSquares)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          IF (in_subset%vertical_levels(idx,block) > 0) THEN
            sumOfSquares    = sumOfSquares + values(idx, block)**2
            ! write(*,*) "number of values=", number_of_values, " value=", values(idx, block), " sum=", sum_value
          ENDIF
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index) &
!ICON_OMP   reduction(+:sumOfSquares)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        sumOfSquares    = sumOfSquares + SUM(values(start_index:end_index, block)**2)
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ENDIF ! (ASSOCIATED(in_subset%vertical_levels))

!     write(0,*) "sumOfSquares=", sumOfSquares
    l2norm = gather_sum(sumOfSquares)
!     write(0,*) "l2norm=", l2norm

    l2norm = SQRT(l2norm)

  END FUNCTION L2Norm_2D_InRange
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the min max mean in a 2D array in a given range subset
  FUNCTION MinMaxMean_2D_InRange(values, in_subset, lacc) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:)
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp) :: minmaxmean(3)
    LOGICAL, INTENT(in), OPTIONAL :: lacc
    LOGICAL :: lzacc

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value, sum_value
    INTEGER :: block, start_index, end_index, number_of_values, idx

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (in_subset%no_of_holes > 0) CALL warning(module_name, "there are holes in the subset")
    ! init the min, max values
    CALL init_min_max(min_value, max_value)
    sum_value = 0._wp
    number_of_values = 0

    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index, idx), reduction(+:number_of_values, sum_value) &
!ICON_OMP reduction(MIN:min_value) reduction(MAX:max_value)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) &
        !$ACC   REDUCTION(+: sum_value, number_of_values) &
        !$ACC   REDUCTION(MAX: max_value) &
        !$ACC   REDUCTION(MIN: min_value) ASYNC(1) IF(lzacc)
        DO idx = start_index, end_index
          IF (in_subset%vertical_levels(idx,block) > 0) THEN
            min_value    = MIN(min_value, values(idx, block))
            max_value    = MAX(max_value, values(idx, block))
            sum_value    = sum_value + values(idx, block)
            number_of_values = number_of_values +  1
            ! write(*,*) "number of values=", number_of_values, " value=", values(idx, block), " sum=", sum_value
          ENDIF
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ELSE ! no in_subset%vertical_levels

#ifdef _OPENACC
        IF (lzacc) CALL finish("MinMaxMean_2D_InRange", "OpenACC version not implemented for global_minmaxmean with no subset")
#endif

!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index, min_in_block, max_in_block) &
!ICON_OMP  reduction(MIN:min_value) reduction(MAX:max_value) reduction(+:sum_value)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        min_in_block = MINVAL(values(start_index:end_index, block))
        max_in_block = MAXVAL(values(start_index:end_index, block))
        min_value    = MIN(min_value, min_in_block)
        max_value    = MAX(max_value, max_in_block)
        sum_value    = sum_value + SUM(values(start_index:end_index, block))
      ENDDO
!ICON_OMP_END_PARALLEL_DO
      ! compute the total number of values
      number_of_values = in_subset%SIZE

    ENDIF ! (ASSOCIATED(in_subset%vertical_levels))

    ! the global min, max, mean, is avaliable only to stdio process
    CALL gather_minmaxmean(min_value, max_value, sum_value, number_of_values, minmaxmean)

  END FUNCTION MinMaxMean_2D_InRange
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the min max mean in a 3D array in a given range subset
  ! The results are over all levels
  FUNCTION MinMaxMean_3D_AllLevels_InRange(values, in_subset, start_level, end_level) result(minmaxmean)
    REAL(wp), INTENT(in) :: values(:,:,:)
    TYPE(t_subset_range), INTENT(in) :: in_subset
    INTEGER, OPTIONAL, INTENT(in) :: start_level, end_level
    REAL(wp) :: minmaxmean(3)

    REAL(wp) :: min_in_block, max_in_block, min_value, max_value, sum_value, global_number_of_values
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical, number_of_values
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':MinMaxMean_3D_AllLevels'

    IF (in_subset%no_of_holes > 0) CALL warning(module_name, "there are holes in the subset")

    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = end_level
    ELSE
      end_vertical = SIZE(values, VerticalDim_Position)
    ENDIF
    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")

    ! init the min, max values
    CALL init_min_max(min_value, max_value)
    sum_value = 0._wp
    number_of_values = 0

    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index, idx, level) &
!ICON_OMP  reduction(MIN:min_value) reduction(MAX:max_value) reduction(+:sum_value, number_of_values)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO idx = start_index, end_index
          !            write(0,*) "end_vertical:", end_vertical," vertical_levels:", in_subset%vertical_levels(idx,block)
          DO level = start_vertical, MIN(end_vertical, in_subset%vertical_levels(idx,block))
            min_value    = MIN(min_value, values(idx, level, block))
            max_value    = MAX(max_value, values(idx, level, block))
            sum_value    = sum_value + values(idx, level, block)
          ENDDO
          number_of_values = number_of_values + MAX(0,    &
            & (MIN(end_vertical, in_subset%vertical_levels(idx,block)) - start_vertical + 1))
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_PARALLEL_DO PRIVATE(block, start_index, end_index, level, min_in_block, max_in_block)  &
!ICON_OMP  reduction(MIN:min_value) reduction(MAX:max_value) reduction(+:sum_value)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        DO level = start_vertical, end_vertical
          min_in_block = MINVAL(values(start_index:end_index, level, block))
          max_in_block = MAXVAL(values(start_index:end_index, level, block))
          min_value    = MIN(min_value, min_in_block)
          max_value    = MAX(max_value, max_in_block)
          sum_value    = sum_value + SUM(values(start_index:end_index, level, block))
        ENDDO
      ENDDO
!ICON_OMP_END_PARALLEL_DO

      IF ((in_subset%end_block - in_subset%start_block) > 1) THEN
        number_of_values = (in_subset%end_block - in_subset%start_block -1) * in_subset%block_size
      ELSE
        number_of_values = 0
      ENDIF

      IF (in_subset%end_block - in_subset%start_block >= 0) THEN
        number_of_values = (number_of_values + in_subset%end_index + &
          & (in_subset%block_size - in_subset%start_index + 1)) * &
          & (end_vertical - start_vertical + 1)
      END IF

    ENDIF

    ! the global min, max, mean, is avaliable only to stdio process
    CALL gather_minmaxmean(min_value, max_value, sum_value, number_of_values, minmaxmean)
    ! write(0,"(a, 1pg26.18, 1pg26.18, 1pg26.18)") "minmaxmean:", minmaxmean(1:3)

  END FUNCTION MinMaxMean_3D_AllLevels_InRange
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  ! Returns the sum in a 3D array in a given indexed subset,
  ! and optional weights with the original 3D index (weights)
  ! or using the the subset index ( subset_indexed_weights ( level, sunset_index)).
  ! The sum is over all levels.
  FUNCTION Sum_3D_AllLevels_3Dweights_InIndexed(values, indexed_subset, start_level, end_level, weights, &
    & subset_indexed_weights) result(total_sum)
    REAL(wp), INTENT(in) :: values(:,:,:)
    TYPE(t_subset_indexed), INTENT(in) :: indexed_subset
    INTEGER,  OPTIONAL, INTENT(in) :: start_level, end_level
    REAL(wp), OPTIONAL, INTENT(in) :: weights(:,:,:)
    REAL(wp), OPTIONAL, INTENT(in) :: subset_indexed_weights(:,:)  ! weights but indexed but the subset index
    ! dim: (vertical_levels, indexed_subset%size)
    REAL(wp) :: total_sum

    REAL(wp) :: sum_value
    INTEGER :: i, block, idx, level,  start_vertical, end_vertical
    INTEGER :: communicator
    !    INTEGER :: idx
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':Sum_3D_AllLevels_3Dweights_InIndexed'


    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = end_level
    ELSE
      end_vertical = SIZE(values, VerticalDim_Position)
    ENDIF
    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")

    ! init the min, max values
    sum_value = 0._wp
    IF (PRESENT(weights)) THEN

      DO i=1, indexed_subset%SIZE
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = start_vertical, end_vertical
          sum_value  = sum_value + values(idx, level, block) * weights(idx, level, block)
        ENDDO
      ENDDO

    ELSEIF (PRESENT(subset_indexed_weights)) THEN

      DO i=1, indexed_subset%SIZE
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = start_vertical, end_vertical
          sum_value  = sum_value + values(idx, level, block) * subset_indexed_weights(level, i)
        ENDDO
      ENDDO

    ELSE

      DO i=1, indexed_subset%SIZE
        block = indexed_subset%block(i)
        idx = indexed_subset%idx(i)
        DO level = start_vertical, end_vertical
          sum_value  = sum_value + values(idx, level, block)
        ENDDO
      ENDDO

    ENDIF

    ! the global min, max is avaliable only to stdio process
    IF (my_process_is_mpi_parallel()) THEN

      communicator = indexed_subset%patch%work_communicator
      ! these are avaliable to all processes
      total_sum = p_sum( sum_value,  comm=communicator)

    ELSE

      total_sum = sum_value

    ENDIF

  END FUNCTION Sum_3D_AllLevels_3Dweights_InIndexed
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Computes the weighted average for each level in a 3D array in a given range subset.
  ! The result is added in the given accumulated_mean array(levels)
  SUBROUTINE AccumulateMean_3D_EachLevel_InRange_2Dweights(values, weights, in_subset, accumulated_mean, start_level, end_level)
    REAL(wp), INTENT(in) :: values(:,:,:) ! in
    REAL(wp), INTENT(in) :: weights(:,:)  ! in
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), INTENT(inout) :: accumulated_mean(:)   ! accumulated mean for each level
    INTEGER, OPTIONAL, INTENT(in) :: start_level, end_level

    REAL(wp) :: mean(SIZE(accumulated_mean))
    INTEGER :: block, level, start_vertical, end_vertical
    INTEGER :: no_of_threads, myThreadNo

    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':AccumulateMean_3D_EachLevel_InRange_2Dweights'

    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = end_level
    ELSE
      end_vertical = SIZE(values, VerticalDim_Position)
    ENDIF
    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")

    CALL LevelHorizontalMean_3D_InRange_2Dweights(values, weights, in_subset, mean, start_vertical, end_vertical)

    ! add average
    DO level = start_vertical, end_vertical
      ! write(0,*) level, ":", total_sum(level), total_weight(level), accumulated_mean(level)
      accumulated_mean(level) = accumulated_mean(level) + mean(level)
    ENDDO

  END SUBROUTINE AccumulateMean_3D_EachLevel_InRange_2Dweights
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the weighted Sum for each level in a 3D array in a given range subset.
  SUBROUTINE LevelHorizontalSum_3D_InRange_2Dweights(values, weights, in_subset, total_sum, start_level, end_level, mean, lopenacc)
    REAL(wp), INTENT(in) :: values(:,:,:) ! in
    REAL(wp), INTENT(in) :: weights(:,:)  ! in
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), INTENT(out) :: total_sum(:)   ! mean for each level
    INTEGER, OPTIONAL, INTENT(in) :: start_level, end_level
    REAL(wp), INTENT(out), OPTIONAL :: mean(:)   ! mean for each level
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc  ! Flag to run on GPU

    ! fixme: this uses a hand-written reduction over threads without
    ! any assertion to be faster/better than the reduction clause of
    ! OpenMP
    REAL(wp), ALLOCATABLE :: sum_value(:,:), sum_weight(:,:), total_weight(:)
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    INTEGER :: allocated_levels, no_of_threads, myThreadNo
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':LevelHorizontalSum_3D_InRange_2Dweights'
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    IF (in_subset%no_of_holes > 0) CALL warning(method_name, "there are holes in the subset")

    no_of_threads = 1
    myThreadNo = 0
#ifdef _OPENMP
    no_of_threads = omp_get_max_threads()
#endif

    allocated_levels = SIZE(total_sum)
    ALLOCATE( sum_value(allocated_levels, 0:no_of_threads-1), &
      & sum_weight(allocated_levels, 0:no_of_threads-1), &
      & total_weight(allocated_levels) )

    !$ACC DATA PRESENT(values, weights, total_sum, mean) &
    !$ACC   CREATE(total_weight) IF(lzopenacc)
    !$ACC DATA CREATE(sum_value, sum_weight) IF(lzopenacc)

    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = end_level
    ELSE
      end_vertical = SIZE(values, VerticalDim_Position)
    ENDIF
    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")

!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
#ifdef _OPENMP
    myThreadNo = omp_get_thread_num()
#endif
!ICON_OMP_SINGLE
#ifdef _OPENMP
    no_of_threads = OMP_GET_NUM_THREADS()
#endif
!ICON_OMP_END_SINGLE NOWAIT
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    !$ACC LOOP GANG VECTOR
    DO level = 1, allocated_levels
      sum_value(level,  myThreadNo) = 0.0_wp
      sum_weight(level, myThreadNo) = 0.0_wp
    END DO
    !$ACC END PARALLEL

    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_DO PRIVATE(block, start_index, end_index, idx, level)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
        !$ACC LOOP SEQ
        DO idx = start_index, end_index
          !$ACC LOOP SEQ
          DO level = start_vertical, MIN(end_vertical, in_subset%vertical_levels(idx,block))
            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & values(idx, level, block) * weights(idx, block)
            sum_weight(level, myThreadNo)  = sum_weight(level, myThreadNo) + weights(idx, block)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_DO PRIVATE(block, start_index, end_index, idx, level)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
        !$ACC LOOP SEQ
        DO idx = start_index, end_index
          ! since we have the same numbder of vertical layers, the weight is the same
          ! for all levels. Compute it only for the first level, and then copy it
          sum_weight(start_vertical, myThreadNo)  = sum_weight(start_vertical, myThreadNo) + weights(idx, block)
          !$ACC LOOP SEQ
          DO level = start_vertical, end_vertical
            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & values(idx, level, block) * weights(idx, block)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

      ! copy the weights to all levels
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
      !$ACC LOOP GANG VECTOR
      DO level = start_vertical+1, end_vertical
         sum_weight(level, myThreadNo)  = sum_weight(start_vertical, myThreadNo)
      ENDDO
      !$ACC END PARALLEL

    ENDIF
!ICON_OMP_END_PARALLEL

    ! gather the total level sum of this process in total_sum(level)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    !$ACC LOOP GANG VECTOR
    DO level = 1, allocated_levels
      total_sum(level)     = 0.0_wp
      total_weight(level) = 0.0_wp
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    !$ACC LOOP SEQ
    DO myThreadNo=0, no_of_threads-1
      !$ACC LOOP GANG VECTOR
      DO level = start_vertical, end_vertical
        ! write(0,*) myThreadNo, level, " sum=", sum_value(level, myThreadNo), sum_weight(level, myThreadNo)
        total_sum(level)    = total_sum(level)    + sum_value(level, myThreadNo)
        total_weight(level) = total_weight(level) + sum_weight(level, myThreadNo)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA
    DEALLOCATE(sum_value, sum_weight)

    ! Collect the value and weight sums (at all procs)
    CALL gather_sums(total_sum, total_weight, lopenacc=lzopenacc)

    IF (PRESENT(mean)) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO level = 1, allocated_levels
        mean(level) = 0.0_wp
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO level = start_vertical, end_vertical
        ! write(0,*) level, ":", total_sum(level), total_weight(level), accumulated_mean(level)
        mean(level) = total_sum(level)/total_weight(level)
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    !$ACC WAIT(1)
    !$ACC END DATA

    DEALLOCATE(total_weight)

  END SUBROUTINE LevelHorizontalSum_3D_InRange_2Dweights
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the weighted average for each level in a 3D array in a given range subset.
  SUBROUTINE LevelHorizontalMean_3D_InRange_2Dweights(values, weights, in_subset, mean, start_level, end_level, lopenacc)
    REAL(wp), INTENT(in) :: values(:,:,:) ! in
    REAL(wp), INTENT(in) :: weights(:,:)  ! in
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), INTENT(out) :: mean(:)   ! mean for each level
    INTEGER, OPTIONAL, INTENT(in) :: start_level, end_level
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc                 ! Flag to run on GPU

    REAL(wp) :: total_sum(SIZE(mean))
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    !$ACC DATA CREATE(total_sum) IF(lzopenacc)

    CALL LevelHorizontalSum_3D_InRange_2Dweights(values=values, weights=weights, in_subset=in_subset, &
      & total_sum=total_sum, start_level=start_level, end_level=end_level, mean=mean, lopenacc=lzopenacc)

    !$ACC END DATA
  END SUBROUTINE LevelHorizontalMean_3D_InRange_2Dweights
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  ! Returns the weighted average for each level in a 3D array in a given range subset.
  REAL(wp) FUNCTION TotalWeightedMean_3D_InRange_3Dweights(values, weights, in_subset, start_level, end_level, lopenacc)

    REAL(wp), INTENT(in) :: values(:,:,:) ! in
    REAL(wp), INTENT(in) :: weights(:,:,:)  ! in
    TYPE(t_subset_range), INTENT(in) :: in_subset
    INTEGER, OPTIONAL, INTENT(in) :: start_level, end_level
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc                 ! Flag to run on GPU

    REAL(wp) :: levelWeights(SIZE(values, VerticalDim_Position)), &
         levelWeightedSum(SIZE(values, VerticalDim_Position))
    REAL(wp) ::  totalWeight, totalSum
    INTEGER :: level, start_vertical, end_vertical
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = end_level
    ELSE
      end_vertical = SIZE(values, VerticalDim_Position)
    ENDIF

    !$ACC DATA CREATE(levelWeights, levelWeightedSum) IF(lzopenacc)

    CALL LevelHorizontalSum_3D_InRange_3Dweights(values=values, weights=weights, in_subset=in_subset, &
      & total_sum=levelWeightedSum, start_level=start_level, end_level=end_level, sumLevelWeights=levelWeights, &
      & lopenacc=lzopenacc )


    totalSum = 0.0_wp
    totalWeight = 0.0_wp
    !$ACC PARALLEL DEFAULT(PRESENT) IF(lzopenacc)
    !$ACC LOOP GANG VECTOR
    DO level = start_vertical, end_vertical
      totalSum    = totalSum    + levelWeightedSum(level)
      totalWeight = totalWeight + levelWeights(level)
    ENDDO
    !$ACC END PARALLEL

    !$ACC END DATA

    TotalWeightedMean_3D_InRange_3Dweights = totalSum / totalWeight

  END FUNCTION TotalWeightedMean_3D_InRange_3Dweights
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the weighted sum for each level in a 3D array in a given range subset.
  SUBROUTINE LevelHorizontalSum_3D_InRange_3Dweights(values, weights, in_subset, total_sum, &
    & start_level, end_level, mean, sumLevelWeights, lopenacc)
    REAL(wp), INTENT(in) :: values(:,:,:) ! in
    REAL(wp), INTENT(in) :: weights(:,:,:)  ! in
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), INTENT(inout) :: total_sum(:)   ! mean for each level
    INTEGER, OPTIONAL, INTENT(in) :: start_level, end_level
    REAL(wp), OPTIONAL, INTENT(out) :: mean(:), sumLevelWeights(:)   ! the sum of the weights in each level
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc                 ! Flag to run on GPU

    ! fixme: this uses a hand-written reduction over threads without
    ! any assertion to be faster/better than the reduction clause of
    ! OpenMP
    REAL(wp), ALLOCATABLE :: sum_value(:,:), sum_weight(:,:), total_weight(:)
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    INTEGER :: allocated_levels, no_of_threads, myThreadNo
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':LevelHorizontalSum_3D_InRange_3Dweights'
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    IF (in_subset%no_of_holes > 0) CALL warning(method_name, "there are holes in the subset")

    no_of_threads = 1
    myThreadNo = 0
#ifdef _OPENMP
    no_of_threads = omp_get_max_threads()
#endif

    allocated_levels = SIZE(total_sum)
    ALLOCATE( sum_value(allocated_levels, 0:no_of_threads-1), &
      & sum_weight(allocated_levels, 0:no_of_threads-1), &
      & total_weight(allocated_levels) )

    !$ACC DATA PRESENT(values, weights, total_sum, mean, sumLevelWeights) &
    !$ACC   CREATE(total_weight) IF(lzopenacc)
    !$ACC DATA CREATE(sum_value, sum_weight) IF(lzopenacc)

    IF (PRESENT(start_level)) THEN
      start_vertical = start_level
    ELSE
      start_vertical = 1
    ENDIF
    IF (PRESENT(end_level)) THEN
      end_vertical = end_level
    ELSE
      end_vertical = SIZE(values, VerticalDim_Position)
    ENDIF
    IF (start_vertical > end_vertical) &
      & CALL finish(method_name, "start_vertical > end_vertical")
    IF ( allocated_levels < end_vertical) &
      & CALL finish(method_name, "allocated_levels < end_vertical")

!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
#ifdef _OPENMP
    myThreadNo = omp_get_thread_num()
#endif
!ICON_OMP_SINGLE
#ifdef _OPENMP
    no_of_threads = OMP_GET_NUM_THREADS()
#endif
!ICON_OMP_END_SINGLE NOWAIT
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    !$ACC LOOP GANG VECTOR
    DO level = 1, allocated_levels
      sum_value(level,  myThreadNo) = 0.0_wp
      sum_weight(level,  myThreadNo) = 0.0_wp
    END DO
    !$ACC END PARALLEL

    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_DO PRIVATE(block, start_index, end_index, idx, level)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
        !$ACC LOOP SEQ
        DO idx = start_index, end_index
          !$ACC LOOP SEQ
          DO level = start_vertical, MIN(end_vertical, in_subset%vertical_levels(idx,block))
            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & values(idx, level, block) * weights(idx, level, block)
            sum_weight(level, myThreadNo)  = sum_weight(level, myThreadNo) + weights(idx, level, block)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_DO PRIVATE(block, start_index, end_index, level)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
        !$ACC LOOP SEQ
        DO idx = start_index, end_index
          ! since we have the same numbder of vertical layers, the weight is the same
          ! for all levels. Compute it only for the first level, and then copy it
          !$ACC LOOP SEQ
          DO level = start_vertical, end_vertical
            sum_value(level, myThreadNo)  = sum_value(level, myThreadNo) + &
              & values(idx, level, block) * weights(idx,level, block)
            sum_weight(level, myThreadNo)  = sum_weight(start_vertical, myThreadNo) + weights(idx, level, block)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

    ENDIF
!ICON_OMP_END_PARALLEL

    ! gather the total level sum of this process in total_sum(level)
    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    !$ACC LOOP GANG VECTOR
    DO level = 1, allocated_levels
      total_sum(level)     = 0.0_wp
      total_weight(level) = 0.0_wp
    END DO
    !$ACC END PARALLEL

    !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    !$ACC LOOP SEQ
    DO myThreadNo=0, no_of_threads-1
      !$ACC LOOP GANG VECTOR
      DO level = start_vertical, end_vertical
        ! write(0,*) myThreadNo, level, " sum=", sum_value(level, myThreadNo), sum_weight(level, myThreadNo)
        total_sum(level)    = total_sum(level)    + sum_value(level, myThreadNo)
        total_weight(level) = total_weight(level) + sum_weight(level, myThreadNo)
      ENDDO
    ENDDO
    !$ACC END PARALLEL
    !$ACC WAIT(1)
    !$ACC END DATA
    DEALLOCATE(sum_value, sum_weight)

    ! Collect the value and weight sums (at all procs)
    CALL gather_sums(total_sum, total_weight, lopenacc=lzopenacc)

    IF (PRESENT(mean)) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO level = 1, allocated_levels
        mean(level) = 0.0_wp
      END DO

      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO level = start_vertical, end_vertical
        mean(level) = total_sum(level)/total_weight(level)
      ENDDO
      !$ACC END PARALLEL
    ENDIF

    IF (PRESENT(sumLevelWeights)) THEN
      !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
      !$ACC LOOP GANG VECTOR
      DO level = 1, allocated_levels
        sumLevelWeights(level) = total_weight(level)
      END DO
      !$ACC END PARALLEL
    ENDIF

    !$ACC WAIT(1)
    !$ACC END DATA

    DEALLOCATE(total_weight)

  END SUBROUTINE LevelHorizontalSum_3D_InRange_3Dweights
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the weighted average for each level in a 3D array in a given range subset.
  SUBROUTINE LevelHorizontalMean_3D_InRange_3Dweights(values, weights, in_subset, mean, start_level, end_level, sumLevelWeights, &
    lopenacc)
    REAL(wp), INTENT(in) :: values(:,:,:) ! in
    REAL(wp), INTENT(in) :: weights(:,:,:)  ! in
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), INTENT(out) :: mean(:)   ! mean for each level
    INTEGER, OPTIONAL :: start_level, end_level
    REAL(wp), OPTIONAL, INTENT(out) :: sumLevelWeights(:)   ! the sum of the weights in each level
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc             ! Flag to run on GPU

    REAL(wp) :: sumLevels(SIZE(mean))
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    !$ACC DATA CREATE(sumLevels) IF(lzopenacc)

    CALL LevelHorizontalSum_3D_InRange_3Dweights(values=values, weights=weights, in_subset=in_subset, &
      & total_sum=sumLevels, start_level=start_level, end_level=end_level, mean=mean, sumLevelWeights=sumLevelWeights, &
      & lopenacc=lzopenacc)

    !$ACC END DATA

  END SUBROUTINE LevelHorizontalMean_3D_InRange_3Dweights
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp)  FUNCTION Sum_2D_InRange(values, in_subset, mean, lopenacc)
    REAL(wp), INTENT(in) :: values(:,:)
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), OPTIONAL, INTENT(out)   :: mean
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc                 ! Flag to run on GPU


    REAL(wp), ALLOCATABLE :: sum_value(:)
    REAL(wp):: total_sum
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    INTEGER :: no_of_threads, myThreadNo, no_of_additions
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':Sum_2D_InRange'
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    IF (in_subset%no_of_holes > 0) CALL warning(module_name, "there are holes in the subset")

    no_of_threads = 1
    myThreadNo = 0
    no_of_additions = 0
#ifdef _OPENMP
    no_of_threads = omp_get_max_threads()
#endif

    ALLOCATE( sum_value(0:no_of_threads-1) )

    !$ACC DATA COPY(values) PRESENT(in_subset, in_subset%vertical_levels) IF(lzopenacc)
    !$ACC DATA CREATE(sum_value) IF(lzopenacc)

!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
#ifdef _OPENMP
    myThreadNo = omp_get_thread_num()
#endif
!ICON_OMP_SINGLE
#ifdef _OPENMP
    no_of_threads = OMP_GET_NUM_THREADS()
#endif
!ICON_OMP_END_SINGLE NOWAIT
    !$ACC KERNELS DEFAULT(PRESENT) IF(lzopenacc)
    sum_value(myThreadNo) = 0.0_wp
    !$ACC END KERNELS
    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_DO PRIVATE(block, start_index, end_index, idx, level)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) IF(lzopenacc)
        !$ACC LOOP SEQ
        DO idx = start_index, end_index
          !$ACC LOOP SEQ
          DO level = 1, MIN(1, in_subset%vertical_levels(idx,block))
            sum_value(myThreadNo) = sum_value(myThreadNo) + values(idx, block)
            no_of_additions = no_of_additions + 1
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_DO PRIVATE(block, start_index, end_index)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) IF(lzopenacc)
        !$ACC LOOP GANG VECTOR
        DO idx = start_index, end_index
          sum_value(myThreadNo)  = sum_value(myThreadNo) + values(idx, block)
          no_of_additions = no_of_additions + 1
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

    ENDIF
!ICON_OMP_END_PARALLEL

    ! gather the total level sum of this process in total_sum(level)
    total_sum     = 0.0_wp
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) &
    !$ACC   REDUCTION(+: total_sum) ASYNC(1) IF(lzopenacc)
    DO myThreadNo=0, no_of_threads-1
      ! write(0,*) myThreadNo, level, " sum=", sum_value(level, myThreadNo)
      total_sum    = total_sum    + sum_value( myThreadNo)
    ENDDO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
    !$ACC END DATA
    DEALLOCATE(sum_value)

    ! Collect the value (at all procs)
    Sum_2D_InRange = gather_sum(total_sum)

    IF (PRESENT(mean)) THEN
      ! Get average
      mean = Sum_2D_InRange / REAL(no_of_additions, KIND=wp)
    ENDIF
    !$ACC END DATA

  END FUNCTION Sum_2D_InRange
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !>
  REAL(wp)  FUNCTION Sum_2D_2Dweights_InRange(values, weights, in_subset, mean, lopenacc)
    REAL(wp), INTENT(in) :: values(:,:)
    REAL(wp), INTENT(in) :: weights(:,:)
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), OPTIONAL, INTENT(out)   :: mean
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc                 ! Flag to run on GPU

    ! fixme: this uses a hand-written reduction over threads without
    ! any assertion to be faster/better than the reduction clause of
    ! OpenMP
    REAL(wp), ALLOCATABLE :: sum_value(:), sum_weight(:)
    REAL(wp):: total_weight, total_sum
    INTEGER :: block, level, start_index, end_index, idx, start_vertical, end_vertical
    INTEGER :: no_of_threads, myThreadNo
    CHARACTER(LEN=*), PARAMETER :: method_name=module_name//':Sum_2D_2Dweights_InRange'
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    IF (in_subset%no_of_holes > 0) CALL warning(module_name, "there are holes in the subset")

    no_of_threads = 1
    myThreadNo = 0
#ifdef _OPENMP
    no_of_threads = omp_get_max_threads()
#endif

    ALLOCATE( sum_value(0:no_of_threads-1), &
      & sum_weight(0:no_of_threads-1) )

    !$ACC DATA PRESENT(values, weights) IF(lzopenacc)
    !$ACC DATA CREATE(sum_value, sum_weight) IF(lzopenacc)

!ICON_OMP_PARALLEL PRIVATE(myThreadNo)
#ifdef _OPENMP
    myThreadNo = omp_get_thread_num()
#endif
!ICON_OMP_SINGLE
#ifdef _OPENMP
    no_of_threads = OMP_GET_NUM_THREADS()
#endif
!ICON_OMP_END_SINGLE NOWAIT
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    sum_value(myThreadNo) = 0.0_wp
    sum_weight(myThreadNo) = 0.0_wp
    !$ACC END KERNELS
    IF (ASSOCIATED(in_subset%vertical_levels)) THEN
!ICON_OMP_DO PRIVATE(block, start_index, end_index, idx, level)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
        !$ACC LOOP SEQ
        DO idx = start_index, end_index
          !$ACC LOOP SEQ
          DO level = 1, MIN(1, in_subset%vertical_levels(idx,block))
            sum_value(myThreadNo)  = sum_value(myThreadNo) + &
              & values(idx, block) * weights(idx, block)
            sum_weight(myThreadNo)  = sum_weight(myThreadNo) + weights(idx, block)
          ENDDO
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

    ELSE ! no in_subset%vertical_levels

!ICON_OMP_DO PRIVATE(block, start_index, end_index)
      DO block = in_subset%start_block, in_subset%end_block
        CALL get_index_range(in_subset, block, start_index, end_index)
        !$ACC PARALLEL DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
        !$ACC LOOP GANG VECTOR
        DO idx = start_index, end_index
          ! since we have the same numbder of vertical layers, the weight is the same
          ! for all levels. Compute it only for the first level, and then copy it
          sum_weight(myThreadNo)  = sum_weight(myThreadNo) + weights(idx, block)
          sum_value(myThreadNo)  = sum_value(myThreadNo) + &
              & values(idx, block) * weights(idx, block)
        ENDDO
        !$ACC END PARALLEL
      ENDDO
!ICON_OMP_END_DO

    ENDIF
!ICON_OMP_END_PARALLEL

    ! gather the total level sum of this process in total_sum(level)
    total_sum    = 0.0_wp
    total_weight = 0.0_wp
    !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) &
    !$ACC   REDUCTION(+: total_sum) REDUCTION(+: total_weight) ASYNC(1) IF(lzopenacc)
    DO myThreadNo=0, no_of_threads-1
      ! write(0,*) myThreadNo, level, " sum=", sum_value(level, myThreadNo), sum_weight(level, myThreadNo)
      total_sum    = total_sum    + sum_value( myThreadNo)
      total_weight = total_weight + sum_weight( myThreadNo)
    ENDDO
    !$ACC END PARALLEL LOOP
    !$ACC WAIT(1)
    !$ACC END DATA
    DEALLOCATE(sum_value, sum_weight)

    ! Collect the value and weight sums (at all procs)
    CALL gather_sums(total_sum, total_weight, lopenacc=lzopenacc)

    IF (PRESENT(mean)) THEN
      ! Get average
      mean = total_sum / total_weight
    ENDIF
    Sum_2D_2Dweights_InRange = total_sum
    !$ACC END DATA
  END FUNCTION Sum_2D_2Dweights_InRange
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  ! Returns the weighted average for each level in a 3D array in a given range subset.
  SUBROUTINE HorizontalMean_2D_InRange_2Dweights(values, weights, in_subset, mean, lopenacc)
    REAL(wp), INTENT(in) :: values(:,:)
    REAL(wp), INTENT(in) :: weights(:,:)
    TYPE(t_subset_range), INTENT(in) :: in_subset
    REAL(wp), INTENT(out) :: mean   ! mean for each level
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc  ! Flag to run on GPU

    REAL(wp) :: WeightedSum
    LOGICAL :: lzopenacc

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    WeightedSum = Sum_2D_2Dweights_InRange(values, weights, in_subset, mean, lopenacc=lzopenacc)

  END SUBROUTINE HorizontalMean_2D_InRange_2Dweights
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  SUBROUTINE init_min_max(min_value, max_value)
    REAL(wp), INTENT(out) :: min_value, max_value

    min_value = 1.e16_wp         ! some large value
    max_value = -1.e16_wp        ! some small value

  END SUBROUTINE init_min_max
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE gather_minmaxmean(min_value, max_value, sum_value, number_of_values, minmaxmean)
    REAL(wp), INTENT(in) :: min_value, max_value, sum_value
    INTEGER, INTENT(in) :: number_of_values
    REAL(wp), INTENT(out) :: minmaxmean(3)

    REAL(wp) :: global_number_of_values
    INTEGER :: communicator

    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      minmaxmean(1) = p_min( min_value,  comm=communicator ) !  mpi_all_reduce
      minmaxmean(2) = p_max( max_value,  comm=communicator )

      ! these are avaliable to all processes
      global_number_of_values = p_sum( REAL(number_of_values,wp),  comm=communicator)
      minmaxmean(3) = p_sum( sum_value,  comm=communicator) / global_number_of_values

    ELSE

      minmaxmean(1) = min_value
      minmaxmean(2) = max_value
      minmaxmean(3) = sum_value / REAL(number_of_values, wp)

    ENDIF
  !  write(0,"(a, 1pg26.18, 1pg26.18, 1pg26.18)") "gather_minmaxmean:", min_value, max_value, sum_value
  !  write(0,*) "number_of_values:",  number_of_values

  END SUBROUTINE gather_minmaxmean
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  REAL(wp) FUNCTION gather_sum(sum_value)
    REAL(wp), INTENT(in) :: sum_value

    INTEGER :: communicator

    IF (my_process_is_mpi_parallel()) THEN
      communicator = get_my_mpi_work_communicator()
      gather_sum = p_sum( sum_value,  comm=communicator)
    ELSE
      gather_sum = sum_value
    ENDIF

  END FUNCTION gather_sum
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE gather_sums_1D(sum_1, sum_2, lopenacc)
    REAL(wp), INTENT(inout) :: sum_1(:), sum_2(:)
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc                 ! Flag to run on GPU

    REAL(wp), ALLOCATABLE :: concat_input_sum(:), concat_output_sum(:)
    INTEGER :: communicator
    INTEGER :: size_of_sum_1, size_of_sum_2, total_size
    LOGICAL :: lzopenacc

    IF (my_process_is_mpi_seq()) RETURN

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    !$ACC UPDATE HOST(sum_1, sum_2) ASYNC(1) IF(lzopenacc)
    !$ACC WAIT(1) IF(lzopenacc)

    size_of_sum_1 = SIZE(sum_1(:))
    size_of_sum_2 = SIZE(sum_2(:))
    total_size   = size_of_sum_1 + size_of_sum_2
    ALLOCATE(concat_input_sum(total_size), concat_output_sum(total_size))

    concat_input_sum(1:size_of_sum_1)            = sum_1(1:size_of_sum_1)
    concat_input_sum(size_of_sum_1+1:total_size) = sum_2(1:size_of_sum_2)
    !write(0,*) "sum_1=", sum_1
    !write(0,*) "sum_2=", sum_2
    !write(0,*) "concat_input_sum=", concat_input_sum

    communicator = get_my_mpi_work_communicator()
    concat_output_sum(:) = p_sum( concat_input_sum,  comm=communicator)

    sum_1(1:size_of_sum_1) = concat_output_sum(1:size_of_sum_1)
    sum_2(1:size_of_sum_2) = concat_output_sum(size_of_sum_1+1:total_size)

    !write(0,*) "concat_output_sum=", concat_output_sum

    DEALLOCATE(concat_input_sum, concat_output_sum)

    !$ACC UPDATE DEVICE(sum_1, sum_2) ASYNC(1) IF(lzopenacc)

  END SUBROUTINE gather_sums_1D
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE gather_sums_0D(sum_1, sum_2, lopenacc)
    REAL(wp), INTENT(inout) :: sum_1, sum_2
    LOGICAL, OPTIONAL, INTENT(in)   :: lopenacc                 ! Flag to run on GPU
    LOGICAL :: lzopenacc

    REAL(wp) :: array_sum_1(1), array_sum_2(1)

    IF (PRESENT(lopenacc)) THEN
      lzopenacc = lopenacc
    ELSE
      lzopenacc = .FALSE.
    ENDIF

    !$ACC DATA CREATE(array_sum_1, array_sum_2) IF(lzopenacc)

    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    array_sum_1(1) = sum_1
    array_sum_2(1) = sum_2
    !$ACC END KERNELS
    CALL gather_sums(array_sum_1, array_sum_2, lopenacc=lzopenacc)
    !$ACC KERNELS DEFAULT(PRESENT) ASYNC(1) IF(lzopenacc)
    sum_1 = array_sum_1(1)
    sum_2 = array_sum_2(1)
    !$ACC END KERNELS

    !$ACC WAIT(1)
    !$ACC END DATA

  END SUBROUTINE gather_sums_0D
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  SUBROUTINE new_statistic_operator(statistic)
    TYPE(t_statistic), INTENT(inout) :: statistic

    statistic%id = new_statistic_no_bars()

  END SUBROUTINE new_statistic_operator
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE delete_statistic_operator(statistic)
    TYPE(t_statistic), INTENT(inout) :: statistic

    CALL delete_statistic(statistic%id)
    statistic%id = -1

  END SUBROUTINE delete_statistic_operator
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new statistics object and returns its id.
  !! The statistic arrays are not allocated.
  INTEGER FUNCTION new_statistic_with_bars(no_of_bars, min_value, max_value, mode)
    INTEGER,  INTENT(in) :: no_of_bars
    REAL(wp), INTENT(in) :: min_value, max_value
    INTEGER,  INTENT(in), OPTIONAL :: mode

    IF (PRESENT(mode)) THEN
      new_statistic_with_bars = new_statistic(mode)
    ELSE
      new_statistic_with_bars = new_statistic(mode)
    ENDIF
    CALL construct_statistic_bars(new_statistic_with_bars, no_of_bars, min_value, max_value)

  END FUNCTION new_statistic_with_bars
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new statistics object and returns its id.
  !! The statistic arrays are not allocated.
  INTEGER FUNCTION new_statistic_no_bars(mode)
    INTEGER,  INTENT(in), OPTIONAL :: mode

    INTEGER :: i

    ! check if the statistic objects have been construced
    IF (no_of_allocated_statistics == 0) THEN
      CALL construct_statistic_objects()
    ENDIF

    IF (max_active_statistics /= active_statistics) THEN
      ! we have a hole (inactive statistic) in the list of max_active_statistics
      DO i = 1, max_active_statistics
        IF (.NOT. statistic_object(i)%is_active) THEN
          new_statistic_no_bars = i
          EXIT
        ENDIF
      ENDDO
    ELSE
      ! add a new statistic object
      IF (max_active_statistics >= no_of_allocated_statistics) THEN
        CALL finish('new_statistic_no_bars', 'exceeded no_of_allocated_statistics');
      ENDIF
      max_active_statistics = max_active_statistics + 1
      new_statistic_no_bars = max_active_statistics
    ENDIF

    active_statistics = active_statistics + 1
    CALL init_statistic_object(new_statistic_no_bars)

    IF (PRESENT(mode)) THEN
      CALL set_mode(new_statistic_no_bars, mode)
    ENDIF


  END FUNCTION new_statistic_no_bars
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new statistics object and returns its id.
  !! The statistic arrays are not allocated.
  INTEGER FUNCTION new_statistic_with_pointers(TARGET,source,mode)
    REAL(wp), TARGET :: TARGET,source
    INTEGER,  INTENT(in), OPTIONAL :: mode
    ! GZ: temporary fix to get the code compiled
    new_statistic_with_pointers = 0
  END FUNCTION new_statistic_with_pointers

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE set_mode(statistic_id, mode)
    INTEGER, INTENT(in) :: statistic_id, mode

    CALL check_active_statistic_id(statistic_id)

    SELECT CASE (mode)

    CASE (add_value, add_max_ratio)
      statistic_object(statistic_id)%mode = mode

    CASE default
      CALL finish('set_mode', 'Unkown mode')

    END SELECT
  END SUBROUTINE set_mode
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deletes the statistic object.
  !! The statistic arrays are deallocated.
  !! The statistic_id is invalid after the delete call.
  !! It maybe associated with another statistic in the future.
  SUBROUTINE delete_statistic(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
    !     CALL deallocate_statistic_object(statistic_id)

    CALL destruct_statistic_bars(statistic_id)
    statistic_object(statistic_id)%is_active = .FALSE.
    active_statistics = active_statistics - 1

  END SUBROUTINE delete_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_data_one_value_real(statistic, val)
    TYPE(t_statistic), INTENT(inout) :: statistic
    REAL(wp), INTENT(in) :: val

    CALL add_statistic_one_value(statistic%id, val)
  END SUBROUTINE add_data_one_value_real
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_data_one_value_int(statistic, val)
    TYPE(t_statistic), INTENT(inout) :: statistic
    INTEGER, INTENT(in) :: val

    CALL add_statistic_one_value(statistic%id, REAL(val,wp))
  END SUBROUTINE add_data_one_value_int
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_statistic_one_value(statistic_id, val)
    INTEGER, INTENT(in) :: statistic_id
    REAL(wp), INTENT(in) :: val

    CALL check_active_statistic_id(statistic_id)

    statistic_object(statistic_id)%sum_of_values  = &
      & statistic_object(statistic_id)%sum_of_values + val
    statistic_object(statistic_id)%max_of_values  = &
      & MAX(statistic_object(statistic_id)%max_of_values, val)
    statistic_object(statistic_id)%min_of_values  = &
      & MIN(statistic_object(statistic_id)%min_of_values, val)

    statistic_object(statistic_id)%no_of_values   = &
      & statistic_object(statistic_id)%no_of_values + 1

  END SUBROUTINE add_statistic_one_value
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE add_statistic_two_values(statistic_id, value1, value2)
    INTEGER, INTENT(in) :: statistic_id
    REAL(wp), INTENT(in) :: value1, value2

    !     REAL(wp) :: max_value, min_value

    !     CALL check_active_statistic_id(statistic_id)

    SELECT CASE (statistic_object(statistic_id)%mode)

    CASE(add_max_ratio)

      CALL add_statistic_one_value(statistic_id, &
        & MAX(value1, value2) / MIN(value1, value2))

    CASE default
      CALL finish('set_mode', 'Unkown mode')

    END SELECT

  END SUBROUTINE add_statistic_two_values
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION min_statistic(statistic)
    TYPE(t_statistic), INTENT(in) :: statistic

    CALL check_active_statistic_id(statistic%id)
    min_statistic = statistic_object(statistic%id)%min_of_values

  END FUNCTION min_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION max_statistic(statistic)
    TYPE(t_statistic), INTENT(in) :: statistic

    CALL check_active_statistic_id(statistic%id)
    max_statistic = statistic_object(statistic%id)%max_of_values

  END FUNCTION max_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION mean_statistic(statistic)
    TYPE(t_statistic), INTENT(in) :: statistic

    CALL check_active_statistic_id(statistic%id)
    mean_statistic = statistic_object(statistic%id)%sum_of_values / &
      & REAL(statistic_object(statistic%id)%no_of_values, wp)

  END FUNCTION mean_statistic
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION min_statistic_of(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
    min_statistic_of = statistic_object(statistic_id)%min_of_values

  END FUNCTION min_statistic_of
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION mean_statistic_of(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
    mean_statistic_of = statistic_object(statistic_id)%sum_of_values / &
      & REAL(statistic_object(statistic_id)%no_of_values, wp)

  END FUNCTION mean_statistic_of
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION max_statistic_of(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    CALL check_active_statistic_id(statistic_id)
    max_statistic_of = statistic_object(statistic_id)%max_of_values

  END FUNCTION max_statistic_of
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Deletes all statistic objects
  !! Note: Should be called only at the stop of a program.
  SUBROUTINE destruct_statistic_objects ()

    INTEGER :: i

    IF (no_of_allocated_statistics == 0) THEN
      CALL warning('destruct_statistic_objects', &
        & 'statistic_objects have not been constructed');
      RETURN
    ENDIF

    DO i=1,max_active_statistics
      IF (statistic_object(i)%is_active) THEN
        CALL delete_statistic(i)
      ENDIF
    ENDDO

    DEALLOCATE(statistic_object)

    no_of_allocated_statistics  = 0
    active_statistics     = 0
    max_active_statistics = 0

  END SUBROUTINE destruct_statistic_objects
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! "Class construction" method.
  !! Called once to initiallize the "class".
  SUBROUTINE construct_statistic_objects()

    INTEGER :: return_status,i

    CHARACTER(*), PARAMETER :: method_name = "construct_statistic_objects"

    IF (no_of_allocated_statistics /= 0) THEN
      ! do some consistency checks, return if everything ok
      IF (no_of_allocated_statistics /= max_no_of_statistic_objects) THEN
        CALL finish(method_name, &
          & 'no_of_allocated_statistics/= max_no_of_statistic_objects');
      ENDIF
      IF (.NOT. ALLOCATED(statistic_object)) THEN
        CALL finish(method_name, '.NOT. ALLOCATED(statistic_object)');
      ENDIF
      IF (SIZE(statistic_object) /= no_of_allocated_statistics) THEN
        CALL finish(method_name, &
          & 'SIZE(statistic_object) /= no_of_allocated_statistics');
      ENDIF

      RETURN
    ENDIF

    IF (ALLOCATED(statistic_object)) THEN
      CALL finish(method_name, 'ALLOCATED(statistic_object)');
    ENDIF

    no_of_allocated_statistics  = max_no_of_statistic_objects
    active_statistics     = 0
    max_active_statistics = 0

    ALLOCATE(statistic_object(no_of_allocated_statistics),stat=return_status)
    IF (return_status /= 0) THEN
      CALL finish(method_name, 'failed to ALLOCATE(statistic_object)');
    ENDIF

    DO i=1,no_of_allocated_statistics
      statistic_object(i)%is_active = .FALSE.
    ENDDO

  END SUBROUTINE construct_statistic_objects
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object
  SUBROUTINE init_statistic_object(statistic_id)
    INTEGER, INTENT(in) :: statistic_id

    statistic_object(statistic_id)%no_of_values   = 0
    statistic_object(statistic_id)%sum_of_values  = 0._wp
    statistic_object(statistic_id)%max_of_values  = &
      & -HUGE(statistic_object(statistic_id)%max_of_values)
    statistic_object(statistic_id)%min_of_values  = &
      & HUGE(statistic_object(statistic_id)%min_of_values)
    statistic_object(statistic_id)%mode     = add_value

    statistic_object(statistic_id)%no_of_bars     = 0
    statistic_object(statistic_id)%min_bars_value = 0._wp
    statistic_object(statistic_id)%max_bars_value = 0._wp
    NULLIFY(statistic_object(statistic_id)%no_of_values_in_bar)

    statistic_object(statistic_id)%is_active = .TRUE.

  END SUBROUTINE init_statistic_object
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object with source and target pointers
  SUBROUTINE init_statistic_object_with_pointers(statistic_id,source,summation,mean,maximum,minimum)
    INTEGER, INTENT(in)            :: statistic_id
    REAL(wp), POINTER :: source
    REAL(wp), POINTER, OPTIONAL :: summation,mean,maximum,minimum

    statistic_object(statistic_id)%no_of_values = 0
    !IF (PRESENT(summation)) THEN
    !  statistic_object(statistic_id)%sum_of_values => summation
    !END IF
    !IF (PRESENT(maximum)) THEN
    !  statistic_object(statistic_id)%max_of_values => maximum
    !END IF
    !IF (PRESENT(minimum)) THEN
    !  statistic_object(statistic_id)%min_of_values => minimum
    !END IF
    statistic_object(statistic_id)%mode = add_value
  END SUBROUTINE init_statistic_object_with_pointers

  !-----------------------------------------------------------------------
  !>
  !! Reset values of the statistic object.
  SUBROUTINE reset_statistic(statistic_id)
    INTEGER, INTENT(in) :: statistic_id
    CALL init_statistic_object(statistic_id)
  END SUBROUTINE reset_statistic

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object
  SUBROUTINE construct_statistic_bars(statistic_id, no_of_bars, min_value, max_value)
    INTEGER,  INTENT(in) :: statistic_id, no_of_bars
    REAL(wp), INTENT(in) :: min_value, max_value

    INTEGER :: return_status
    CHARACTER(*), PARAMETER :: method_name = "construct_statistic_bars"

    CALL check_active_statistic_id(statistic_id)
    IF (no_of_bars < 2) &
      & CALL finish(method_name,'no_of_bars < 2');
    IF (min_value >= max_value) &
      & CALL finish(method_name,'min_value >= max_value');

    statistic_object(statistic_id)%no_of_bars     = no_of_bars
    statistic_object(statistic_id)%min_bars_value = min_value
    statistic_object(statistic_id)%max_bars_value = max_value

    ALLOCATE(statistic_object(statistic_id)%no_of_values_in_bar(no_of_bars),stat=return_status)
    IF (return_status /= 0) THEN
      CALL finish(method_name, 'failed to ALLOCATE(statistic_object)');
    ENDIF

  END SUBROUTINE construct_statistic_bars
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Initializes a new statistic object
  SUBROUTINE destruct_statistic_bars(statistic_id)
    INTEGER,  INTENT(in) :: statistic_id

    !     INTEGER :: return_status
    !     CHARACTER(*), PARAMETER :: method_name = "destruct_statistic_bars"

    CALL check_active_statistic_id(statistic_id)

    IF (ASSOCIATED(statistic_object(statistic_id)%no_of_values_in_bar)) THEN
      DEALLOCATE(statistic_object(statistic_id)%no_of_values_in_bar)
    ENDIF

    statistic_object(statistic_id)%no_of_bars     = 0
    NULLIFY(statistic_object(statistic_id)%no_of_values_in_bar)

  END SUBROUTINE destruct_statistic_bars
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Checks if a the statistic object associated with the statistic_id is active.
  !! If the statistic object is not active the program stops with an error.
  SUBROUTINE check_active_statistic_id (statistic_id)
    INTEGER :: statistic_id

    IF (statistic_id > max_active_statistics) THEN
      CALL finish('get_statistic', 'statistic_id > max_active_statistics');
    ENDIF
    IF (.NOT. statistic_object(statistic_id)%is_active) THEN
      CALL finish('get_statistic', 'statistic is not active');
    ENDIF
  END SUBROUTINE check_active_statistic_id
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE add_fields_3d(sum_field,field,subset,levels,has_missvals, missval)
    REAL(wp),INTENT(inout)          :: sum_field(:,:,:)
    REAL(wp),INTENT(in)             :: field(:,:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    INTEGER,INTENT(in),OPTIONAL :: levels
    LOGICAL, INTENT(IN), OPTIONAL :: has_missvals
    REAL(wp), INTENT(IN), OPTIONAL :: missval

    INTEGER :: idx,block,level,start_index,end_index

    INTEGER :: mylevels
    LOGICAL :: my_force_level, my_has_missvals
    REAL(wp) :: my_miss

    my_has_missvals = .FALSE.
    my_miss = 0.0_wp
    IF (PRESENT(has_missvals)) my_has_missvals = has_missvals
    IF (PRESENT(missval)) my_miss = missval

    ! use constant levels
    mylevels                       = SIZE(sum_field, VerticalDim_Position)
    IF (PRESENT(levels))  mylevels = levels
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
    DO block = subset%start_block, subset%end_block
      CALL get_index_range(subset, block, start_index, end_index)
      !$ACC PARALLEL PRESENT(sum_field, field) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR COLLAPSE(2)
      DO level = 1, mylevels
        DO idx = start_index, end_index
          sum_field(idx,level,block) = MERGE(my_miss, &
                                        & sum_field(idx,level,block) + field(idx,level,block), &
                                        & my_has_missvals .AND. (field(idx,level,block) == my_miss))
        END DO
      END DO
      !$ACC END PARALLEL
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE add_fields_3d

  SUBROUTINE add_fields_2d(sum_field,field,subset,has_missvals, missval)
    REAL(wp),INTENT(inout)          :: sum_field(:,:)
    REAL(wp),INTENT(in)             :: field(:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    LOGICAL, INTENT(IN), OPTIONAL :: has_missvals
    REAL(wp), INTENT(IN), OPTIONAL :: missval

    INTEGER :: jb,jc,start_index,end_index
    LOGICAL :: my_has_missvals
    REAL(wp) :: my_miss

    my_has_missvals = .FALSE.
    my_miss = 0.0_wp
    IF (PRESENT(has_missvals)) my_has_missvals = has_missvals
    IF (PRESENT(missval)) my_miss = missval

!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc) SCHEDULE(dynamic)
    DO jb = subset%start_block, subset%end_block
      CALL get_index_range(subset, jb, start_index, end_index)
      !$ACC PARALLEL PRESENT(sum_field, field) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = start_index, end_index
        sum_field(jc,jb) = MERGE(my_miss, sum_field(jc,jb) + field(jc,jb), my_has_missvals .AND. (field(jc,jb) == my_miss))
      END DO
      !$ACC END PARALLEL
    END DO
!ICON_OMP_END_PARALLEL_DO
  END SUBROUTINE add_fields_2d

  SUBROUTINE add_fields_3d_sp(sum_field,field,subset,levels,has_missvals, missval)
    REAL(wp),INTENT(inout)          :: sum_field(:,:,:)
    REAL(sp),INTENT(in)             :: field(:,:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    INTEGER,INTENT(in),OPTIONAL :: levels
    LOGICAL, INTENT(IN), OPTIONAL :: has_missvals
    REAL(wp), INTENT(IN), OPTIONAL :: missval

    INTEGER :: idx,block,level,start_index,end_index

    INTEGER :: mylevels
    LOGICAL :: my_force_level, my_has_missvals
    REAL(wp) :: my_miss

    my_has_missvals = .FALSE.
    my_miss = 0.0_wp
    IF (PRESENT(has_missvals)) my_has_missvals = has_missvals
    IF (PRESENT(missval)) my_miss = missval

      ! use constant levels
      mylevels                       = SIZE(sum_field, VerticalDim_Position)
      IF (PRESENT(levels))  mylevels = levels
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        !$ACC PARALLEL PRESENT(sum_field, field) ASYNC(1) IF(i_am_accel_node)
        !$ACC LOOP GANG VECTOR COLLAPSE(2)
        DO level = 1, mylevels
          DO idx = start_index, end_index
            sum_field(idx,level,block) = MERGE(my_miss, &
                                             & sum_field(idx,level,block) + REAL(field(idx,level,block),wp), &
                                             & my_has_missvals .AND. (field(idx,level,block) == REAL(my_miss,sp)))
          END DO
        END DO
        !$ACC END PARALLEL
      END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE add_fields_3d_sp

  SUBROUTINE add_fields_2d_sp(sum_field,field,subset,has_missvals, missval)
    REAL(wp),INTENT(inout)          :: sum_field(:,:)
    REAL(sp),INTENT(in)             :: field(:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    LOGICAL, INTENT(IN), OPTIONAL :: has_missvals
    REAL(wp), INTENT(IN), OPTIONAL :: missval

    INTEGER :: jb,jc,start_index,end_index
    LOGICAL :: my_has_missvals
    REAL(wp) :: my_miss

    my_has_missvals = .FALSE.
    my_miss = 0.0_wp
    IF (PRESENT(has_missvals)) my_has_missvals = has_missvals
    IF (PRESENT(missval)) my_miss = missval

!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, jc) SCHEDULE(dynamic)
    DO jb = subset%start_block, subset%end_block
      CALL get_index_range(subset, jb, start_index, end_index)
      !$ACC PARALLEL PRESENT(sum_field, field) ASYNC(1) IF(i_am_accel_node)
      !$ACC LOOP GANG VECTOR
      DO jc = start_index, end_index
        sum_field(jc,jb) = MERGE(my_miss, &
                               & sum_field(jc,jb) + REAL(field(jc,jb),wp), &
                               & my_has_missvals .AND. (field(jc,jb) == REAL(my_miss,sp)))
      END DO
      !$ACC END PARALLEL
    END DO
!ICON_OMP_END_PARALLEL_DO
  END SUBROUTINE add_fields_2d_sp

  !-----------------------------------------------------------------------
  SUBROUTINE add_verticallyIntegrated_field(vint_field_acc,field_3D,subset,height,levels)
    REAL(wp),INTENT(inout)          :: vint_field_acc(:,:)
    REAL(wp),INTENT(in)             :: field_3D(:,:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    REAL(wp),INTENT(in)             :: height(:,:,:)
    INTEGER,INTENT(in),OPTIONAL :: levels

    INTEGER :: idx,block,level,start_index,end_index

    INTEGER :: mylevels
    LOGICAL :: my_force_level

    IF (ASSOCIATED(subset%vertical_levels) .AND. .NOT. PRESENT(levels)) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = 1, subset%vertical_levels(idx,block)
            vint_field_acc(idx,block) = vint_field_acc(idx,block) + field_3D(idx,level,block) * height(idx,level,block)
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSE
      ! use constant levels
      mylevels   = SIZE(field_3D, VerticalDim_Position)
      IF (PRESENT(levels)) mylevels = levels
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level)
!SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = 1, mylevels
            vint_field_acc(idx,block) = vint_field_acc(idx,block) + field_3D(idx,level,block) * height(idx,level,block)
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ENDIF
  END SUBROUTINE add_verticallyIntegrated_field

  !-----------------------------------------------------------------------
  SUBROUTINE verticallyIntegrated_field_weighted(vint_field,field_3D,subset,height,levels)
    REAL(wp),INTENT(inout)          :: vint_field(:,:)
    REAL(wp),INTENT(in)             :: field_3D(:,:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    REAL(wp),INTENT(in)             :: height(:,:,:)
    INTEGER,INTENT(in),OPTIONAL :: levels

    INTEGER :: idx,block,level,start_index,end_index

    INTEGER :: mylevels
    LOGICAL :: my_force_level

    IF (ASSOCIATED(subset%vertical_levels) .AND. .NOT. PRESENT(levels)) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        DO idx = start_index, end_index
          vint_field(idx,block) = 0.0_wp
          DO level = 1, subset%vertical_levels(idx,block)
            vint_field(idx,block) = vint_field(idx,block) + field_3D(idx,level,block) * height(idx,level,block)
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSE
      ! use constant levels
      mylevels   = SIZE(field_3D, VerticalDim_Position)
      IF (PRESENT(levels)) mylevels = levels
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        DO idx = start_index, end_index
           vint_field(idx,block) = 0.0_wp
           DO level = 1, mylevels
            vint_field(idx,block) = vint_field(idx,block) + field_3D(idx,level,block) * height(idx,level,block)
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ENDIF
  END SUBROUTINE verticallyIntegrated_field_weighted
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE verticallyIntegrated_field_notWeighted(vint_field,field_3D,subset,levels, lacc)
    REAL(wp),INTENT(inout)          :: vint_field(:,:)
    REAL(wp),INTENT(in)             :: field_3D(:,:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    INTEGER,INTENT(in),OPTIONAL :: levels
    LOGICAL, INTENT(IN), OPTIONAL :: lacc

    INTEGER :: idx,block,level,start_index,end_index

    INTEGER :: mylevels
    LOGICAL :: my_force_level
    LOGICAL  :: lzacc

    CALL set_acc_host_or_device(lzacc, lacc)

    IF (ASSOCIATED(subset%vertical_levels) .AND. .NOT. PRESENT(levels)) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO idx = start_index, end_index
          vint_field(idx,block) = 0.0_wp
          !$ACC LOOP SEQ
          DO level = 1, subset%vertical_levels(idx,block)
            vint_field(idx,block) = vint_field(idx,block) + field_3D(idx,level,block)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSE
      ! use constant levels
      mylevels   = SIZE(field_3D, VerticalDim_Position)
      IF (PRESENT(levels)) mylevels = levels
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) IF(lzacc)
        DO idx = start_index, end_index
           vint_field(idx,block) = 0.0_wp
          !$ACC LOOP SEQ
           DO level = 1, mylevels
            vint_field(idx,block) = vint_field(idx,block) + field_3D(idx,level,block)
          END DO
        END DO
        !$ACC END PARALLEL LOOP
      END DO
!ICON_OMP_END_PARALLEL_DO

    ENDIF
  END SUBROUTINE verticallyIntegrated_field_notWeighted
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE add_verticalSum_field(vint_field_acc,field_3D,subset,levels)
    REAL(wp),INTENT(inout)          :: vint_field_acc(:,:)
    REAL(wp),INTENT(in)             :: field_3D(:,:,:)
    TYPE(t_subset_range),INTENT(in) :: subset
    INTEGER,INTENT(in),OPTIONAL :: levels

    INTEGER :: idx,block,level,start_index,end_index

    INTEGER :: mylevels
    LOGICAL :: my_force_level

    IF (ASSOCIATED(subset%vertical_levels) .AND. .NOT. PRESENT(levels)) THEN
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = 1, subset%vertical_levels(idx,block)
            vint_field_acc(idx,block) = vint_field_acc(idx,block) + field_3D(idx,level,block)
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ELSE
      ! use constant levels
      mylevels   = SIZE(field_3D, VerticalDim_Position)
      IF (PRESENT(levels)) mylevels = levels
!ICON_OMP_PARALLEL_DO PRIVATE(start_index, end_index, idx, level) SCHEDULE(dynamic)
      DO block = subset%start_block, subset%end_block
        CALL get_index_range(subset, block, start_index, end_index)
        DO idx = start_index, end_index
          DO level = 1, mylevels
            vint_field_acc(idx,block) = vint_field_acc(idx,block) + field_3D(idx,level,block)
          END DO
        END DO
      END DO
!ICON_OMP_END_PARALLEL_DO

    ENDIF
  END SUBROUTINE add_verticalSum_field
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Computes updated time average
  !!
  !! Computes updated time average for a particular field
  !!
  !! @Literature
  !! Based on proposal found in
  !! Jochen Froehlich, 2006:Large Eddy Simulation turbulenter Stroemungen, Teubner,
  !! page 273
  !!
  FUNCTION time_avg (psi_avg_old, psi_inst, wgt)  RESULT (psi_avg_new)

    REAL(wp), INTENT(IN) :: psi_avg_old       !< time average at t(n-1)
    REAL(wp), INTENT(IN) :: psi_inst          !< instantaneous value
    REAL(wp), INTENT(IN) :: wgt               !< weight (=dt/sim_time)

    ! Result
    REAL(wp) :: psi_avg_new                   !< updated time average

#ifdef _OPENACC
    !$ACC ROUTINE SEQ
#endif

    !--------------------------------------------------------------------

    ! compute updated time average
    psi_avg_new = (1._wp - wgt)*psi_avg_old + wgt*psi_inst

  END FUNCTION time_avg

END MODULE mo_statistics
