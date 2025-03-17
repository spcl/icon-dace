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

MODULE mo_dist_dir

  USE mo_exception, ONLY: finish
  USE mo_mpi, ONLY: p_alltoall, p_alltoallv

  IMPLICIT NONE

  TYPE t_dist_dir
    INTEGER, ALLOCATABLE :: owner(:)
    INTEGER :: global_size, local_start_index
    INTEGER :: comm, comm_size, comm_rank
  END TYPE t_dist_dir

  PUBLIC :: dist_dir_setup
  PUBLIC :: dist_dir_get_owners
  PUBLIC :: t_dist_dir

  INTERFACE dist_dir_get_owners
    MODULE PROCEDURE dist_dir_get_owners_all
    MODULE PROCEDURE dist_dir_get_owners_subset
  END INTERFACE dist_dir_get_owners

CONTAINS

  !> Function for setting up a distributed directory
  !> @param[in] owned_indices array containing all indices owned by the
  !>                          local process (indices have to be in the range
  !>                          1 : global_size)
  !> @param[in] global_size   highest global index
  !> @param[in] comm          communicator that contains all processes that
  !>                          are part of the distributed directory
  !> @param[in] comm_rank     rank of process in the provided communicator
  !> @param[in] comm_size     number of processes in the provided
  !>                          communicator
  !> @returns distributed directory
  !> @remarks this function is collective for all processes in comm
  SUBROUTINE dist_dir_setup(dist_dir, owned_indices, global_size, comm, &
    &                       comm_rank, comm_size)

    TYPE(t_dist_dir), INTENT(OUT) :: dist_dir
    INTEGER, INTENT(IN) :: owned_indices(:)
    INTEGER, INTENT(IN) :: global_size
    INTEGER, INTENT(IN) :: comm
    INTEGER, INTENT(IN) :: comm_rank
    INTEGER, INTENT(IN) :: comm_size

    INTEGER :: num_indices_per_process
    INTEGER, POINTER :: recv_buffer(:)
    INTEGER :: dummy(comm_size)
    INTEGER :: num_recv_indices_per_process(comm_size)
    INTEGER :: i, j, n

    num_indices_per_process = (global_size + comm_size - 1) / comm_size
    dist_dir%global_size = global_size
    dist_dir%local_start_index = comm_rank * num_indices_per_process + 1
    dist_dir%comm = comm
    dist_dir%comm_rank = comm_rank
    dist_dir%comm_size = comm_size

    ALLOCATE(dist_dir%owner(num_indices_per_process))

    CALL distribute_indices(owned_indices, recv_buffer, dummy, &
      &                     num_recv_indices_per_process, &
      &                     num_indices_per_process, comm, comm_size)

    dist_dir%owner(:) = -1
    n = 0
    DO i = 1, comm_size
      DO j = 1, num_recv_indices_per_process(i)
        n = n + 1
        dist_dir%owner(recv_buffer(n)) = i - 1
      END DO
    END DO

    DEALLOCATE(recv_buffer)

  END SUBROUTINE dist_dir_setup

  !> gets for each provided global index the rank of the process it
  !> belongs to
  !> @param[in] dist_dir distributed directory
  !> @param[in] indices  indices for which the owner are to be retrieved
  !> @return returns the owners of the provided indices
  !> @remark this routine is collective for all processes that belong to
  !>         dist_dir
  !> @remark all provided global indices need to be valid
  FUNCTION dist_dir_get_owners_all(dist_dir, indices) RESULT(owners)

    TYPE(t_dist_dir), INTENT(IN) :: dist_dir
    INTEGER, INTENT(IN) :: indices(:)
    INTEGER :: owners(SIZE(indices))

    INTEGER :: num_indices_per_process
    INTEGER :: num_send_indices_per_process(dist_dir%comm_size)
    INTEGER :: num_recv_indices_per_process(dist_dir%comm_size)
    INTEGER, POINTER :: recv_buffer(:)
    INTEGER, ALLOCATABLE :: send_buffer(:)
    INTEGER :: send_displ(dist_dir%comm_size), recv_displ(dist_dir%comm_size+1)

    INTEGER :: i, n

    IF (ANY(indices(:) < 1 .OR. indices(:) > dist_dir%global_size)) &
      CALL finish("get_owner_ranks", "invalid global index")

    recv_buffer => null()

    num_indices_per_process = (dist_dir%global_size + dist_dir%comm_size - 1) &
      &                       / dist_dir%comm_size

    CALL distribute_indices(indices, recv_buffer, &
      &                     num_recv_indices_per_process, &
      &                     num_send_indices_per_process, &
      &                     num_indices_per_process, dist_dir%comm, &
      &                     dist_dir%comm_size)

    ALLOCATE(send_buffer(SIZE(recv_buffer(:))))

    send_buffer(:) = dist_dir%owner(recv_buffer(:))

    DEALLOCATE(recv_buffer)
    ALLOCATE(recv_buffer(SIZE(indices(:))))

    send_displ(1) = 0
    recv_displ(1) = 0
    DO i = 2, dist_dir%comm_size
      send_displ(i) = send_displ(i-1) + num_send_indices_per_process(i-1)
      recv_displ(i) = recv_displ(i-1) + num_recv_indices_per_process(i-1)
    END DO
    recv_displ(dist_dir%comm_size+1) = recv_displ(dist_dir%comm_size) + &
      num_recv_indices_per_process(dist_dir%comm_size)

    CALL p_alltoallv(send_buffer, num_send_indices_per_process, send_displ, &
      &              recv_buffer, num_recv_indices_per_process, recv_displ, &
      &              dist_dir%comm)

    DO i = SIZE(indices(:)), 1, -1
      n = 2 + (indices(i)-1) / num_indices_per_process
      owners(i) = recv_buffer(recv_displ(n))
      recv_displ(n) = recv_displ(n) - 1
    END DO

    DEALLOCATE(send_buffer, recv_buffer)

  END FUNCTION dist_dir_get_owners_all

  !> gets for each provided global index of indices where mask is also true
  !! the rank of the process it belongs to, -1 otherwise
  !!
  !! @param[in] dist_dir distributed directory
  !! @param[in] indices  indices for which the owner are to be retrieved
  !! @param[in] mask only fetch owner(i) where mask(i) is .TRUE.
  !> @return returns the owners of the provided indices
  !! @remark this routine is collective for all processes that belong to
  !!         dist_dir
  !! @remark all provided global indices need to be valid
  FUNCTION dist_dir_get_owners_subset(dist_dir, indices, mask) &
       RESULT(owners)
    TYPE(t_dist_dir), INTENT(IN) :: dist_dir
    LOGICAL, INTENT(in) :: mask(:)
    INTEGER :: owners(1:SIZE(mask))
    INTEGER, INTENT(in) :: indices(:)

    INTEGER, ALLOCATABLE :: temp(:)



    INTEGER :: m, n


    m = COUNT(mask)
    n = SIZE(mask)

    IF (SIZE(indices) /= n .OR. n /= SIZE(owners) .OR. m > n) &
         CALL finish("get_owner_subset_from_dist_dir", "invalid arguments")

    ALLOCATE(temp(m+1))
    owners(1:m) = PACK(indices(:), mask(:))
    temp(1:m) = dist_dir_get_owners(dist_dir, owners(1:m))
    ! xlf generates a BPT trap in the unpack code when compiled with -qcheck

    owners(1:n) = UNPACK(temp(1:m), mask(1:n), -1)







  END FUNCTION dist_dir_get_owners_subset

  SUBROUTINE distribute_indices(indices, recv_buffer, &
    &                           num_send_indices_per_process, &
    &                           num_recv_indices_per_process, &
    &                           num_indices_per_process, comm, comm_size)

    INTEGER, INTENT(IN) :: indices(:)
    INTEGER, INTENT(OUT), POINTER :: recv_buffer(:)
    INTEGER, INTENT(IN) :: comm, comm_size
    INTEGER, INTENT(OUT) :: num_send_indices_per_process(comm_size)
    INTEGER, INTENT(OUT) :: num_recv_indices_per_process(comm_size)
    INTEGER, INTENT(IN) :: num_indices_per_process

    INTEGER :: send_displ(comm_size+1), recv_displ(comm_size)
    INTEGER, ALLOCATABLE :: send_buffer(:)
    INTEGER :: i, j

    num_send_indices_per_process(:) = 0
    DO i = 1, SIZE(indices(:))
      j = 1 + (indices(i)-1) / num_indices_per_process
      num_send_indices_per_process(j) = &
        num_send_indices_per_process(j) + 1
    END DO

    call p_alltoall(num_send_indices_per_process(:), &
      &             num_recv_indices_per_process(:), comm)

    ALLOCATE(send_buffer(SUM(num_send_indices_per_process(:))), &
      &      recv_buffer(SUM(num_recv_indices_per_process(:))))

    send_displ(1:2) = 0
    recv_displ(1) = 0
    DO i = 2, comm_size
      send_displ(i+1) = send_displ(i) + num_send_indices_per_process(i-1)
      recv_displ(i)   = recv_displ(i-1) + num_recv_indices_per_process(i-1)
    END DO

    DO i = 1, SIZE(indices(:))
      j = 2 + (indices(i)-1) / num_indices_per_process
      send_displ(j) = send_displ(j) + 1
      send_buffer(send_displ(j)) = MOD(indices(i)-1, num_indices_per_process) + 1
    END DO

    CALL p_alltoallv(send_buffer, num_send_indices_per_process, send_displ, &
      &              recv_buffer, num_recv_indices_per_process, recv_displ, &
      &              comm)

    DEALLOCATE(send_buffer)

  END SUBROUTINE distribute_indices

END MODULE mo_dist_dir
