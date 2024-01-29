!> Contains some functions to interface with the parallel infrastructure
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2024, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_jsb_parallel
#ifndef __NO_JSBACH__

  USE mo_jsb_parallel_iface
#ifdef _OPENMP
  USE omp_lib, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
#endif

  IMPLICIT NONE
  PUBLIC

CONTAINS

  LOGICAL FUNCTION Is_enabled_openmp()

#ifdef _OPENMP
    Is_enabled_openmp = .TRUE.
#else
    Is_enabled_openmp = .FALSE.
#endif

  END FUNCTION Is_enabled_openmp

  LOGICAL FUNCTION Is_omp_inside_serial()

    IF (Is_enabled_openmp()) THEN
#ifdef _OPENMP
      Is_omp_inside_serial = omp_get_num_threads() == 1
#else
      Is_omp_inside_serial = .TRUE.
#endif
    ELSE
      Is_omp_inside_serial = .TRUE.
    END IF

  END FUNCTION Is_omp_inside_serial

  INTEGER FUNCTION Get_omp_no_of_threads()

    IF (Is_enabled_openmp()) THEN
#ifdef _OPENMP
      Get_omp_no_of_threads = omp_get_max_threads()
#else
      Get_omp_no_of_threads = 1
#endif
    ELSE
      Get_omp_no_of_threads = 1
    END IF

  END FUNCTION Get_omp_no_of_threads

  INTEGER FUNCTION Get_omp_thread()

    IF (Is_enabled_openmp()) THEN
#ifdef _OPENMP
      Get_omp_thread = omp_get_thread_num() + 1
#else
      Get_omp_thread = 1
#endif
    ELSE
      Get_omp_thread = 1
    END IF

  END FUNCTION Get_omp_thread

#endif
END MODULE mo_jsb_parallel