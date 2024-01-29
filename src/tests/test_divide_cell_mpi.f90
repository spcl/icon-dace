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

PROGRAM test_divide_cell_mpi

#ifndef NOMPI
  USE mo_decomposition_tools, ONLY: t_cell_info, &
       divide_cells_by_location, sort_cell_info_by_cell_number, &
       OPERATOR(/=)
  USE mo_mpi, ONLY: start_mpi, stop_mpi, abort_mpi, &
       process_mpi_all_comm, mpi_in_place, mpi_success, p_bool, mpi_lor, &
       split_global_mpi_communicator, p_int, p_pe_work, p_n_work, &
       p_barrier
  USE mo_divide_cells_by_location_mpi, ONLY: divide_cells_by_location_mpi, &
       init_divide_cells_by_location_mpi
  USE mpi
#endif
  USE mo_io_units, ONLY: nerr

  IMPLICIT NONE

#ifdef NOMPI
  WRITE (nerr, '(a)') 'MPI test skipped in no-MPI configuration.'
#else

  LOGICAL, PARAMETER :: debug = .FALSE.
  INTEGER :: ierror
  INTEGER :: comm
  INTEGER :: comm_rank, comm_size
  CHARACTER(*), PARAMETER :: method_name = 'test_divide_cell_mpi'
  INTEGER, ALLOCATABLE :: rand_init_data(:)


  CALL start_mpi('test_divide_cell_mpi')
  CALL split_global_mpi_communicator(1, 1)

  comm = process_mpi_all_comm
  CALL mpi_comm_size(comm, comm_size, ierror)
  IF (ierror /= mpi_success) THEN
    WRITE (nerr,'(a,a)') method_name, ' mpi_comm_size failed.'
    WRITE (nerr,'(a,i4)') ' Error =  ', ierror
    CALL abort_mpi
  END IF
  CALL mpi_comm_rank(comm, comm_rank, ierror)
  IF (ierror /= mpi_success) THEN
    WRITE (nerr,'(a,a)') method_name, ' mpi_comm_rank failed.'
    WRITE (nerr,'(a,i4)') ' Error =  ', ierror
    CALL abort_mpi
  END IF

  CALL init_random_state
  CALL init_divide_cells_by_location_mpi

  CALL randomized_comparison()
  CALL stop_mpi

CONTAINS
  ELEMENTAL FUNCTION cell_desc2str(c) RESULT(s)
    TYPE(t_cell_info), INTENT(in) :: c
    INTEGER, PARAMETER :: idig = RANGE(1) + 1 + 1
    CHARACTER(4 * idig + (4 - 1) * 2) :: s
    WRITE (s, '(3(i0,", "),i0)') c%lat, c%lon, c%cell_number, c%owner
  END FUNCTION cell_desc2str

  SUBROUTINE comparison_run(cell_desc, nparts, firstpart)
    TYPE(t_cell_info), INTENT(in) :: cell_desc(:)
    INTEGER, INTENT(in) :: nparts, firstpart

    TYPE(t_cell_info), ALLOCATABLE :: cell_desc_par(:), cell_desc_ser(:)

    INTEGER :: ncells_g
    INTEGER :: my_input_start, my_input_end, my_part_start, my_part_end
    INTEGER :: i, owner, ierror
    LOGICAL, ALLOCATABLE :: seen(:)
    CHARACTER(*), PARAMETER :: method_name &
         = 'test_divide_cell_mpi::comparison_run'

    ncells_g = SIZE(cell_desc)
    my_input_start = (comm_rank * ncells_g) / comm_size + 1
    my_input_end = ((comm_rank + 1) * ncells_g) / comm_size
    my_part_start = (comm_rank * nparts) / comm_size
    my_part_end = ((comm_rank + 1) * nparts) / comm_size - 1

    ! the MAX prevents zero-size allocations problems with compilers like xlf
    ! in -qzerosize mode
    ALLOCATE(cell_desc_par(my_input_start:&
         MAX(my_input_end, my_input_start+1)), &
         cell_desc_ser(1:ncells_g))

    cell_desc_ser(:) = cell_desc(:)
    IF (my_input_end >= my_input_start) &
      cell_desc_par(my_input_start:my_input_end) &
           = cell_desc(my_input_start:my_input_end)

    CALL divide_cells_by_location_mpi(ncells_g, &
         my_input_end - my_input_start + 1, &
         cell_desc_par(:), nparts, comm)

    CALL divide_cells_by_location(ncells_g, cell_desc_ser, &
         firstpart, nparts - 1 + firstpart)
    CALL p_barrier
    CALL sort_cell_info_by_cell_number(cell_desc_ser, ncells_g)

    ALLOCATE(seen(firstpart:nparts + firstpart))
    seen = .FALSE.
    seen(firstpart:nparts+firstpart-1) = .FALSE.
    DO i = my_input_start, my_input_end
      owner = cell_desc_par(i)%owner
      IF (owner < firstpart .OR. owner >= nparts + firstpart) THEN
        WRITE (nerr, '(a,i0,2a)') ' cell ', i, &
             'found invalid owner ASSIGNMENT!', &
             TRIM(cell_desc2str(cell_desc_par(i)))
        CALL test_abort
      END IF
      seen(owner) = .TRUE.
    END DO
    ! fixme: use something scalable here, perhaps based on mpi_reduce_scatter
    CALL mpi_allreduce(mpi_in_place, seen, ncells_g, p_bool, mpi_lor, comm, ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_allreduce failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL test_abort
    END IF
    DO i = my_part_start, my_part_end
      IF (.NOT. seen(i)) THEN
        WRITE (nerr, '(a,i0,a)') 'no cells assigned to part ', i, '!'
        CALL test_abort
      END IF
    END DO

    DO i = my_input_start, my_input_end
      IF (cell_desc_ser(i) /= cell_desc_par(i)) THEN
        WRITE (nerr, '(a,i0,4a)') 'reference mismatch at position ', i, &
             ': ', TRIM(cell_desc2str(cell_desc_ser(i))), '; ', &
             TRIM(cell_desc2str(cell_desc_par(i)))
      END IF
    END DO

    IF (debug) THEN
      PRINT *, ''
      DO i = 1, ncells_g
        PRINT *, TRIM(cell_desc2str(cell_desc_ser(i)))
      END DO
      IF (comm_rank == 0) THEN
        PRINT *, ''
        DO i = 1, ncells_g
          PRINT *, TRIM(cell_desc2str(cell_desc_ser(i)))
        END DO
      END IF
    END IF
    DEALLOCATE(seen)
    DEALLOCATE(cell_desc_par)
    DEALLOCATE(cell_desc_ser)
  END SUBROUTINE comparison_run

  SUBROUTINE init_random_state
    INTEGER, PARAMETER :: n_date_values=8
    INTEGER :: i, m, n
    CALL RANDOM_SEED(size=n)
    m = MAX(n, n_date_values)
    ALLOCATE(rand_init_data(m))
    rand_init_data(:) = 0
    CALL DATE_AND_TIME(values=rand_init_data(1:n_date_values))
    DO i = n_date_values, n
      rand_init_data(i)  = rand_init_data(i - n_date_values + 1)
    END DO
    DO i = 1, n
      rand_init_data(i) = rand_init_data(i) + rand_init_data(n - i + 1)
    END DO
    CALL mpi_bcast(rand_init_data, n_date_values, p_int, 0, comm, ierror)
    IF (ierror /= mpi_success) THEN
      WRITE (nerr,'(a,a)') method_name, ' mpi_bcast failed.'
      WRITE (nerr,'(a,i4)') ' Error =  ', ierror
      CALL test_abort
    END IF
    IF (debug .AND. p_pe_work == 0) &
         WRITE (0, '(a,i0,(7(",",i0)))') 'rand_init_data=', &
         rand_init_data(1:n)
    CALL RANDOM_SEED(put=rand_init_data(1:n))
  END SUBROUTINE init_random_state

  SUBROUTINE simple_comparison()
    INTEGER, PARAMETER :: ncells = 10
    INTEGER, PARAMETER :: nparts = 3, firstpart = 0
    TYPE(t_cell_info), PARAMETER :: cell_desc_ref(ncells) &
         = (/ t_cell_info(-5000, -1000, 1, HUGE(1)), &
         &    t_cell_info(5000, -1000, 2, -HUGE(1)), &
         &    t_cell_info(-5000, 1000, 3, HUGE(1)), &
         &    t_cell_info(5000, 1000, 4, -HUGE(1)), &
         &    t_cell_info(-5000, 0, 5, HUGE(1)), &
         &    t_cell_info(-5000, -2000, 6, HUGE(1)), &
         &    t_cell_info(5000, -2000, 7, -HUGE(1)), &
         &    t_cell_info(-5000, 2000, 8, HUGE(1)), &
         &    t_cell_info(5000, 2000, 9, -HUGE(1)), &
         &    t_cell_info(-5000, 0, 10, HUGE(1)) /)
    CALL comparison_run(cell_desc_ref, nparts, firstpart)
  END SUBROUTINE simple_comparison

  SUBROUTINE randomized_comparison()
    USE mo_math_utilities, ONLY: fxp_lat, fxp_lon
    USE mo_math_constants, ONLY: pi
    USE iso_c_binding, ONLY: c_int
    TYPE(t_cell_info), ALLOCATABLE :: cells(:)
    REAL :: r(2)
    INTEGER :: i, j, n, nparts, firstpart
    INTERFACE
      SUBROUTINE usleep(usecs) BIND(c, name='usleep')
        IMPORT :: c_int
        INTEGER(c_int), VALUE :: usecs
      END SUBROUTINE usleep
    END INTERFACE

    CALL RANDOM_NUMBER(r)
    n = CEILING(r(1) * 100 * comm_size)
    nparts = CEILING(r(2) * 10)
    IF (nparts > n) n = nparts + NINT(r(1))
    ALLOCATE(cells(n))
    firstpart = 0
    DO i = 1, n
      CALL RANDOM_NUMBER(r)
      cells(i)%lon = fxp_lon(2 * pi * (r(1) - 0.5))
      cells(i)%lat = fxp_lat(pi * (r(2) - 0.5))
      cells(i)%cell_number = i
      CALL RANDOM_NUMBER(r)
      ! check if this cell should duplicate another's lat/lon or both values
      IF (r(1) < 0.5**4) THEN
        j = INT(r(2) * (i - 1)) + 1
        IF (r(1) < 0.5**5) THEN
          ! copy lon
          cells(i)%lon = cells(j)%lon
        ELSE IF (r(1) >= 0.5**5) THEN
          ! copy lat
          cells(i)%lat = cells(j)%lat
        END IF
      END IF
    END DO
    IF (debug) THEN
      DO j = 1, n, 2
        DO i = 0, p_n_work - 1
          IF (p_pe_work == i) THEN
            WRITE (0, '(i0,3(a,", "))') p_pe_work, ': ', &
                 cell_desc2str(cells(j:MIN(j+1, n)))
          END IF
          CALL p_barrier
        END DO
        CALL usleep(50000_c_int)
      END DO
    END IF
    CALL comparison_run(cells, nparts, firstpart)
  END SUBROUTINE randomized_comparison

  SUBROUTINE test_abort
    INTEGER, PARAMETER :: rand_seed_unit = 10
    IF (p_pe_work == 0) THEN
      OPEN(unit=rand_seed_unit, file='rand_seed.txt', status='replace', &
           form='formatted')
      WRITE (unit=rand_seed_unit, fmt='(9(i0,", "),i0)') &
           rand_init_data
      CLOSE(rand_seed_unit)
    END IF
    CALL abort_mpi
  END SUBROUTINE test_abort
#endif
END PROGRAM test_divide_cell_mpi
