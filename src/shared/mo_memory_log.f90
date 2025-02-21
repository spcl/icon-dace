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

MODULE mo_memory_log
  use mo_mpi, only: get_my_mpi_all_id, get_my_mpi_all_comm_size
  use mo_util_rusage
  use mo_util_sysinfo

  IMPLICIT NONE

  LOGICAL, SAVE :: memory_log_active = .false.
  INTEGER :: memHandle
  INTEGER :: memory_log_id

  ! public routines
  PUBLIC :: memory_log_initialize
  PUBLIC :: memory_log_terminate
  PUBLIC :: memory_log_add
  PUBLIC :: memory_log_print_maxrss

  contains

    subroutine memory_log_initialize(log_all_ranks)
      LOGICAL, OPTIONAL :: log_all_ranks

      CHARACTER(len=10) ::  tag
      LOGICAL           :: my_log_all_ranks

      my_log_all_ranks = .FALSE.

      IF (PRESENT(log_all_ranks)) my_log_all_ranks = log_all_ranks

      IF (my_log_all_ranks) THEN
        memory_log_active = .true.
      ELSE
        ! catch the root process
        IF (0 .eq. get_my_mpi_all_id()) THEN
          memory_log_active = .true.
          ! the highest rank
        ELSE IF ( (get_my_mpi_all_comm_size()-1) .eq. get_my_mpi_all_id() ) THEN
          memory_log_active = .true.
          ! and something in-between
        ELSE IF ( (get_my_mpi_all_comm_size()/2) .eq. get_my_mpi_all_id() ) THEN
          memory_log_active = .true.
        ENDIF
      ENDIF

      IF (.NOT. memory_log_active) RETURN




      memory_log_id = 0


      write(tag,'(i8.8)') memory_log_id
      memHandle = add_rss_list('memUsage',tag=tag)
    end subroutine memory_log_initialize

    subroutine memory_log_terminate
      IF (.NOT. memory_log_active ) RETURN
      CALL close_rss_lists()
    end subroutine memory_log_terminate

    subroutine memory_log_add
      IF (.NOT. memory_log_active ) RETURN
      CALL add_rss_usage(memHandle)
    end subroutine memory_log_add

    subroutine memory_log_print_maxrss
      INTEGER :: maxrss

      CALL util_get_maxrss(maxrss)
      PRINT  *, "PE #", memory_log_id,": MAXRSS (MiB) = ", maxrss
    end subroutine memory_log_print_maxrss

end module
