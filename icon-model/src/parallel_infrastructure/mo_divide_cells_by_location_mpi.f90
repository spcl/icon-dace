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

MODULE mo_divide_cells_by_location_mpi
#ifndef NOMPI
  USE mo_kind, ONLY: i8, wp
  USE mo_mpi, ONLY: mpi_op_null, mpi_datatype_null, mpi_success, abort_mpi, &
       mpi_in_place, p_int, p_int_i8, mpi_address_kind, mpi_sum, mpi_2integer, &
       mpi_integer
  USE mpi, ONLY: MPI_Op_create
#ifndef NO_MPI_CHOICE_ARG
  USE mpi, ONLY: MPI_Allreduce, MPI_Get_address
#endif
#ifdef SLOW_MPI_USER_REDUCTION
  USE mpi, ONLY: mpi_min, mpi_max
#endif
  USE mo_math_utilities, ONLY: flp_lon, flp_lat
  USE mo_decomposition_tools, ONLY: t_cell_info
  USE mo_io_units, ONLY: nerr

  IMPLICIT NONE
  PRIVATE

  INTEGER :: divide_cells_reduce_op = mpi_op_null, &
       cell_num_bounds_reduce_op = mpi_op_null
  INTEGER :: divide_cells_reduce_dt = mpi_datatype_null
  TYPE t_divide_cells_agg
    SEQUENCE
    INTEGER :: max_lat, max_lon, min_lat, min_lon
    INTEGER(i8) :: lat_sum, lon_sum
  END TYPE t_divide_cells_agg

  PUBLIC :: init_divide_cells_by_location_mpi
  PUBLIC :: divide_cells_by_location_mpi

CONTAINS

  SUBROUTINE init_divide_cells_by_location_mpi
#ifndef SLOW_MPI_USER_REDUCTION
    IF (divide_cells_reduce_dt == mpi_datatype_null) CALL create_dc_op_and_dt
#endif
  END SUBROUTINE init_divide_cells_by_location_mpi

#ifndef SLOW_MPI_USER_REDUCTION
  SUBROUTINE divide_cells_agg(a, b, n, dt)
    INTEGER, INTENT(in) :: n, dt
    TYPE(t_divide_cells_agg), INTENT(in) :: a(n)
    TYPE(t_divide_cells_agg), INTENT(inout) :: b(n)

    INTEGER :: i
    DO i = 1, n
      b(i)%min_lat = MIN(b(i)%min_lat, a(i)%min_lat)
      b(i)%max_lat = MAX(b(i)%max_lat, a(i)%max_lat)
      b(i)%min_lon = MIN(b(i)%min_lon, a(i)%min_lon)
      b(i)%max_lon = MAX(b(i)%max_lon, a(i)%max_lon)
      b(i)%lat_sum = b(i)%lat_sum + a(i)%lat_sum
      b(i)%lon_sum = b(i)%lon_sum + a(i)%lon_sum
    END DO
  END SUBROUTINE divide_cells_agg

  SUBROUTINE cell_num_bounds_agg(a, b, n, dt)
    INTEGER, INTENT(in) :: n, dt
    INTEGER, INTENT(in) :: a(2, n)
    INTEGER, INTENT(inout) :: b(2, n)

    INTEGER :: i
    DO i = 1, n
      b(1, i) = MIN(b(1, i), a(1, i))
      b(2, i) = MAX(b(2, i), a(2, i))
    END DO
  END SUBROUTINE cell_num_bounds_agg

  SUBROUTINE create_dc_op_and_dt
    CHARACTER(len=*), PARAMETER :: method_name = 'create_dc_op'
    INTEGER, PARAMETER :: nblk = 2
    INTEGER, PARAMETER :: dt_bl(nblk) = (/ 4, 2 /)
    INTEGER :: dt_comp_dt(nblk)
    INTEGER(mpi_address_kind) :: dt_displ(nblk)
    TYPE(t_divide_cells_agg) :: dummy
    INTEGER :: i, ierror
    CALL mpi_op_create(divide_cells_agg, .TRUE., divide_cells_reduce_op, &
         ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_create_op failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF

    CALL mpi_op_create(cell_num_bounds_agg, .TRUE., cell_num_bounds_reduce_op, &
         ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_create_op failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF

    dt_comp_dt(1) = p_int
    dt_comp_dt(2) = p_int_i8
    CALL mpi_get_address(dummy, dt_displ(1), ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_get_address failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF
    CALL mpi_get_address(dummy%lat_sum, dt_displ(2), ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_get_address failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF
    DO i = nblk, 1, -1
      dt_displ(i) = dt_displ(i) - dt_displ(1)
    END DO
    CALL mpi_type_create_struct(nblk, dt_bl, dt_displ, dt_comp_dt, &
         divide_cells_reduce_dt, ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_type_create_struct failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF
    CALL mpi_type_commit(divide_cells_reduce_dt, ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_type_commit failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF
  END SUBROUTINE create_dc_op_and_dt
#endif

  !-------------------------------------------------------------------------
  !>
  !! Divide along lat/lon axes in parallel
  !!
  SUBROUTINE divide_cells_by_location_mpi(n_cells_g, n_cells, cell_desc, &
       npart, comm)
    INTEGER, INTENT(in) :: n_cells_g, n_cells, npart, comm

    CHARACTER(len=*), PARAMETER :: method_name = 'divide_cells_by_location_par'
    TYPE(t_cell_info), INTENT(inout) :: cell_desc(:)

    INTEGER :: depth, selmask, path
    INTEGER :: target_part_size, part_ofs
    INTEGER, PARAMETER :: dm_lat = 0, dm_lon = 2, dm_lat_num = 1, dm_lon_num = 3
    TYPE t_pivot
      INTEGER :: lat_or_lon, cell_number, div_method
    END TYPE t_pivot
    TYPE stack
      INTEGER :: nparts, ncells
      TYPE(t_pivot) :: pivot
    END TYPE stack
    TYPE(t_pivot) :: divider
    INTEGER :: nsubpart, i, ierror
    TYPE(stack) :: part_size_stack(0:BIT_SIZE(1))
    INTEGER :: cell_num_bounds(2)
#ifdef SLOW_MPI_IN_PLACE
    INTEGER :: cell_num_bounds_temp(2)
#endif

    cell_num_bounds(1) = HUGE(cell_num_bounds(1))
    cell_num_bounds(2) = -HUGE(cell_num_bounds(2))
    ! establish minimal/maximal cell numbers to bound later interval nesting
    DO i = 1, SIZE(cell_desc)
      cell_num_bounds(1) = MIN(cell_num_bounds(1), cell_desc(i)%cell_number)
      cell_num_bounds(2) = MAX(cell_num_bounds(2), cell_desc(i)%cell_number)
    END DO
#ifndef SLOW_MPI_USER_REDUCTION
    CALL mpi_allreduce(mpi_in_place, cell_num_bounds, 1, MPI_2INTEGER, &
         cell_num_bounds_reduce_op, comm, ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF
#else
#  ifdef SLOW_MPI_IN_PLACE
    CALL mpi_allreduce(cell_num_bounds(1), cell_num_bounds_temp(1), &
         1, mpi_integer, mpi_min, comm, ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF
    CALL mpi_allreduce(cell_num_bounds(2), cell_num_bounds_temp(2), &
         1, mpi_integer, mpi_max, comm, ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL abort_mpi
    END IF
    cell_num_bounds = cell_num_bounds_temp
#  else
    "mpi_in_place without mpi user reductions unimplemented"
#  endif
#endif

    ! set owner to 0 to reflect intermediate assignment to root set
    ! at first the corresponding path will be temporarily stored in %owner
    ! to denote which cells to process at each step
    cell_desc(:)%owner = 0
    ! construct root case
    ! depth is 0, but set to 1 for first descent
    depth = 1
    ! at the start the path records no decision
    path = 0
    ! choose bit which identifies cells on current depth
    selmask = IBSET(0, BIT_SIZE(1) - depth)
    ! we are finished if the recorded path is all ones and no further
    ! descent is necessary

    part_size_stack(0)%nparts = npart
    part_size_stack(0)%ncells = n_cells_g

    nsubpart = npart
    target_part_size = INT(DBLE(n_cells_g)*DBLE(npart/2)/DBLE(npart))

    preorder_part_walk : DO WHILE (depth > 0)
      ! in case no sub partitioning has to be performed, the current
      ! path can be backtracked, i.e. ascent until finding a node
      ! where the right child hasn't been processed
      DO WHILE (nsubpart == 1)
        depth = depth - 1
        IF (depth == 0) EXIT preorder_part_walk
        DO WHILE (IAND(path, IBSET(0, BIT_SIZE(1) - depth)) /= 0)
          path = IBCLR(path, BIT_SIZE(1) - depth)
          depth = depth - 1
          IF (depth == 0) EXIT preorder_part_walk
        END DO
        path = IBSET(path, BIT_SIZE(1) - depth)
        nsubpart = part_size_stack(depth - 1)%nparts
        target_part_size = part_size_stack(depth - 1)%ncells
        target_part_size = target_part_size &
             - INT(DBLE(target_part_size) * DBLE(nsubpart/2) / DBLE(nsubpart))
        nsubpart = nsubpart - nsubpart/2
        part_size_stack(depth)%ncells = target_part_size
        part_size_stack(depth)%nparts = nsubpart
        target_part_size = INT(DBLE(target_part_size) &
             * DBLE(nsubpart/2) / DBLE(nsubpart))
        depth = depth + 1
      END DO
      ! find pivot such that part_size_a elements are <= pivot and
      ! part_size_b elements are > pivot
      ! this also sets the comparison method if comparison by lat/lon is
      ! insufficient to generate the desired bisection
      divider = find_pivot()
      ! first note appended path for each cell
      DO i = 1, n_cells
        IF (cell_desc(i)%owner == path) THEN
          cell_desc(i)%owner = IOR(path, &
               ISHFT(pivot_cmp(cell_desc(i), divider), BIT_SIZE(1) - depth))
        END IF
      END DO
      ! given the divider, this function performs an iterative
      ! preorder traversal of the partitions

      ! set up for processing of left branch of descent
      part_size_stack(depth)%ncells = target_part_size
      ! how many partitions to generate for "left" branch of conceptual
      ! tree walk
      nsubpart = nsubpart / 2
      part_size_stack(depth)%nparts = nsubpart
      target_part_size = INT(DBLE(target_part_size) &
           &                 * DBLE(nsubpart / 2) / DBLE(nsubpart))
      depth = depth + 1

      ! processing of the right branch of descent happens implicitly through
      ! the ascend operation above
    END DO preorder_part_walk

    ! rewrite owner field to partition number from stored path
    DO i = 1, n_cells
      nsubpart = npart
      path = cell_desc(i)%owner
      part_ofs = 0
      DO depth = BIT_SIZE(1) - 1, 0, -1
        selmask = IBSET(0, depth)
        IF (IAND(path, selmask) /= 0) THEN
          part_ofs = part_ofs + nsubpart / 2
          nsubpart = nsubpart - nsubpart / 2
        ELSE
          nsubpart = nsubpart / 2
        END IF
      END DO
      cell_desc(i)%owner = part_ofs
    END DO

  CONTAINS
    FUNCTION find_pivot() RESULT(pivot)
      TYPE(t_pivot) :: pivot

      CHARACTER(len=*), PARAMETER :: method_name &
           = 'divide_cells_by_location_par::find_pivot'
      TYPE(t_divide_cells_agg) :: agg
      INTEGER :: i, ierror
      INTEGER(i8) :: nom1, denom1, nom2, denom2
      INTEGER, PARAMETER :: max_pivot_le = 64
      ! candidates for median
      INTEGER :: pivot_le(max_pivot_le, 2)
      INTEGER :: npivot_le, ncells2divide
      INTEGER :: pivot_guess_min, pivot_guess_med, pivot_guess_max
      REAL(wp) :: min_lat, max_lat, min_lon, max_lon, avglat, scalexp
#ifdef SLOW_MPI_IN_PLACE
      TYPE(t_divide_cells_agg) :: agg_temp
      INTEGER :: pivot_le_temp(max_pivot_le)
#endif
      agg%min_lat = HUGE(agg%min_lat)
      agg%max_lat = -HUGE(agg%max_lat)
      agg%min_lon = HUGE(agg%min_lon)
      agg%max_lon = -HUGE(agg%max_lon)
      agg%lat_sum = 0_i8
      agg%lon_sum = 0_i8

      DO i = 1, n_cells
        IF (cell_desc(i)%owner == path) THEN
          agg%min_lat = MIN(agg%min_lat, cell_desc(i)%lat)
          agg%max_lat = MAX(agg%max_lat, cell_desc(i)%lat)
          agg%min_lon = MIN(agg%min_lon, cell_desc(i)%lon)
          agg%max_lon = MAX(agg%max_lon, cell_desc(i)%lon)
          agg%lat_sum = agg%lat_sum + cell_desc(i)%lat
          agg%lon_sum = agg%lon_sum + cell_desc(i)%lon
        END IF
      END DO

#ifndef SLOW_MPI_USER_REDUCTION
      CALL mpi_allreduce(mpi_in_place, agg, 1, divide_cells_reduce_dt, &
           divide_cells_reduce_op, comm, ierror)
      IF (ierror /= mpi_success) THEN
        WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', ierror
        CALL abort_mpi
      END IF
#else
#  ifdef SLOW_MPI_IN_PLACE
      CALL mpi_allreduce(agg%min_lat, agg_temp%min_lat, 2, mpi_integer, &
           mpi_min, comm, ierror)
      IF (ierror /= mpi_success) THEN
        WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', ierror
        CALL abort_mpi
      END IF
      CALL mpi_allreduce(agg%max_lat, agg_temp%max_lat, 2, mpi_integer, &
           mpi_max, comm, ierror)
      IF (ierror /= mpi_success) THEN
        WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', ierror
        CALL abort_mpi
      END IF
      CALL mpi_allreduce(agg%lat_sum, agg_temp%lat_sum, 2, p_int_i8, &
           mpi_sum, comm, ierror)
      IF (ierror /= mpi_success) THEN
        WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
        WRITE (nerr,'(a,i4)') ' Error =  ', ierror
        CALL abort_mpi
      END IF
      agg = agg_temp
#  else
     "mpi_in_place without mpi user reductions unimplemented"
#  endif
#endif

      min_lon = flp_lon(agg%min_lon)
      min_lat = flp_lat(agg%min_lat)
      max_lon = flp_lon(agg%max_lon)
      max_lat = flp_lat(agg%max_lat)
      ! average latitude in patch
      ncells2divide = part_size_stack(depth-1)%ncells
      avglat = flp_lat(INT(agg%lat_sum/INT(ncells2divide,i8)))
      ! account somehow for convergence of meridians - this formula is just empiric
      scalexp = 1._wp - MAX(0._wp,ABS(max_lat)-1._wp,ABS(min_lat)-1._wp)
      min_lon = min_lon*(COS(avglat))**scalexp
      max_lon = max_lon*(COS(avglat))**scalexp

      IF (max_lat - min_lat >= max_lon - min_lon) THEN
        pivot%div_method = dm_lat
        pivot_guess_min = agg%min_lat
        pivot_guess_med = INT(agg%lat_sum / INT(ncells2divide, i8))
        pivot_guess_max = agg%max_lat
      ELSE
        pivot%div_method = dm_lon
        pivot_guess_min = agg%min_lon
        pivot_guess_med = INT(agg%lon_sum / INT(ncells2divide, i8))
        pivot_guess_max = agg%max_lon
      END IF

      npivot_le = MAX(MIN(max_pivot_le, ncells2divide), 3)

      find_pivot_le: DO WHILE (.TRUE.)
        pivot_le(npivot_le/2 + 1, 1) = pivot_guess_med
        pivot_le(npivot_le/2 + 1, 2) = count_le(pivot, pivot_guess_med, pivot%div_method, path)

        nom1   = INT(pivot_guess_med - pivot_guess_min, i8)
        denom1 = INT(npivot_le/2, i8)
        nom2   = INT(pivot_guess_max - pivot_guess_med, i8)
        denom2 = INT(npivot_le - npivot_le/2 - 1, i8)
!$OMP PARALLEL
!$OMP DO
        DO i = 1, npivot_le/2
          pivot_le(i, 1) = pivot_guess_min + INT(INT(i - 1, i8) * nom1 / denom1)
          pivot_le(i, 2) = count_le(pivot, pivot_le(i, 1), pivot%div_method, path)
        END DO
!$OMP END DO NOWAIT
!$OMP DO
        DO i = npivot_le/2 + 2, npivot_le
          pivot_le(i, 1) = pivot_guess_med + INT(INT(i - npivot_le/2 - 1, i8) * nom2 / denom2)
          pivot_le(i, 2) = count_le(pivot, pivot_le(i, 1), pivot%div_method, path)
        END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

        CALL mpi_allreduce(mpi_in_place, pivot_le(1:npivot_le, 2), &
             npivot_le, mpi_integer, mpi_sum, comm, ierror)
        IF (ierror /= mpi_success) THEN
          WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
          WRITE (nerr,'(a,i4)') ' Error =  ', ierror
          CALL abort_mpi
        END IF
        DO i = 1, npivot_le
          IF (pivot_le(i, 2) == target_part_size) THEN
            IF (IAND(pivot%div_method, 1) == 0) THEN
              pivot%lat_or_lon = pivot_le(i, 1)
            ELSE
              pivot%cell_number = pivot_le(i, 1)
            END IF
            RETURN
          END IF
          ! once we are well beyond the target there is no need to search further
          IF (pivot_le(i, 2) > target_part_size) EXIT
        END DO
        ! are more than half of coords close to min?
        ! then cell_number needs to be used for disambiguation
        IF (pivot_le(1, 2) > target_part_size) THEN
          IF (.NOT. pivot%div_method == dm_lat &
               .AND. .NOT. pivot%div_method == dm_lon) THEN
            CALL abort_mpi
          END IF
          pivot%lat_or_lon = pivot_le(1, 1)
          pivot%div_method = pivot%div_method + 1
          pivot_guess_min = cell_num_bounds(1)
          pivot_guess_max = cell_num_bounds(2)
          pivot_guess_med = (cell_num_bounds(1) + cell_num_bounds(2)) / 2
        ELSE
          ! pivot_le(i, 2) > target_part_size
          pivot_guess_min = pivot_le(i - 1, 1)
          pivot_guess_med = pivot_le(i - 1, 1) &
               + (pivot_le(i, 1) - pivot_le(i - 1, 1))/2
          pivot_guess_max = pivot_le(i, 1)
          ! are pivot_guess_min and _max so close together, that no new
          ! pivot range results? (i.e. both pivot candidates have
          ! already been tested)
          IF (pivot_guess_min == pivot_guess_med) THEN
            ! this must not happen for cell_number based discrimation since
            ! no two cells have the same cell_number and thus
            IF (.NOT. pivot%div_method == dm_lat &
                 .AND. .NOT. pivot%div_method == dm_lon) THEN
              CALL abort_mpi
            END IF
            pivot%lat_or_lon = pivot_guess_max
            pivot%div_method = pivot%div_method + 1
            pivot_guess_min = cell_num_bounds(1)
            pivot_guess_max = cell_num_bounds(2)
            pivot_guess_med = (cell_num_bounds(1) + cell_num_bounds(2)) / 2
          END IF
        END IF
      END DO find_pivot_le
    END FUNCTION find_pivot

    ELEMENTAL FUNCTION count_le(pivot, pivot_val, method, ipath)
      TYPE(t_pivot), INTENT(in) :: pivot
      INTEGER, INTENT(in) :: pivot_val, method, ipath
      INTEGER :: count_le
      SELECT CASE (method)
      CASE (dm_lat)
        count_le = COUNT(cell_desc%owner == ipath &
             &           .AND. cell_desc%lat .LE. pivot_val)
      CASE (dm_lon)
        count_le = COUNT(cell_desc%owner == ipath &
             &           .AND. cell_desc%lon .LE. pivot_val)
      CASE (dm_lat_num)
        count_le = COUNT(cell_desc%owner == ipath &
             &           .AND. (cell_desc%lat < pivot%lat_or_lon &
             &                  .OR. (cell_desc%lat == pivot%lat_or_lon &
             &                        .AND. cell_desc%cell_number &
             &                              .LE. pivot_val)))
      CASE (dm_lon_num)
        count_le = COUNT(cell_desc%owner == ipath &
             &           .AND. (cell_desc%lon < pivot%lat_or_lon &
             &                  .OR. (cell_desc%lon == pivot%lat_or_lon &
             &                        .AND. cell_desc%cell_number &
             &                              .LE. pivot_val)))
      END SELECT
    END FUNCTION count_le

    ELEMENTAL FUNCTION pivot_cmp(cell_desc, pivot) RESULT(cmp)
      TYPE(t_cell_info), INTENT(in) :: cell_desc
      TYPE(t_pivot), INTENT(in) :: pivot
      INTEGER :: cmp
      SELECT CASE (pivot%div_method)
      CASE (dm_lat)
        cmp = MERGE(1, 0, cell_desc%lat > pivot%lat_or_lon)
      CASE (dm_lon)
        cmp = MERGE(1, 0, cell_desc%lon > pivot%lat_or_lon)
      CASE (dm_lat_num)
        cmp = MERGE(1, 0, cell_desc%lat > pivot%lat_or_lon &
             &            .OR. (cell_desc%lat == pivot%lat_or_lon &
             &                  .AND. cell_desc%cell_number &
             &                              > pivot%cell_number))
      CASE (dm_lon_num)
        cmp = MERGE(1, 0, cell_desc%lon > pivot%lat_or_lon &
             &            .OR. (cell_desc%lon == pivot%lat_or_lon &
             &                  .AND. cell_desc%cell_number &
             &                              > pivot%cell_number))
      END SELECT
    END FUNCTION pivot_cmp

  END SUBROUTINE divide_cells_by_location_mpi

#endif
END MODULE mo_divide_cells_by_location_mpi
