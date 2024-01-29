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

MODULE mo_util_timer
  USE ISO_C_BINDING, ONLY: c_double, c_int, c_ptr

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: util_cputime
  PUBLIC :: util_walltime
  PUBLIC :: util_gettimeofday
  PUBLIC :: util_init_real_time
  PUBLIC :: util_get_real_time_size
  PUBLIC :: util_read_real_time
  PUBLIC :: util_diff_real_time

  INTERFACE
    FUNCTION util_cputime(user_time, system_time) BIND(C) RESULT(f_result)
      IMPORT c_double, c_int
      REAL(c_double), INTENT(OUT) :: user_time, system_time
      INTEGER(c_int) :: f_result
    END FUNCTION util_cputime

    FUNCTION util_walltime() BIND(C) RESULT(f_result)
      IMPORT c_double
      REAL(c_double) :: f_result
    END FUNCTION util_walltime

    FUNCTION util_gettimeofday() BIND(C) RESULT(f_result)
      IMPORT c_double
      REAL(c_double) :: f_result
    END FUNCTION util_gettimeofday

    SUBROUTINE util_init_real_time() BIND(C)
    END SUBROUTINE util_init_real_time

    SUBROUTINE util_get_real_time_size(rt_size) BIND(C)
      IMPORT c_int
      INTEGER(c_int), INTENT(OUT) :: rt_size
    END SUBROUTINE util_get_real_time_size

    SUBROUTINE util_read_real_time(it) BIND(C)
      IMPORT c_ptr
      TYPE(c_ptr), VALUE :: it
    END SUBROUTINE util_read_real_time

    SUBROUTINE util_diff_real_time(it1, it2, t) BIND(C)
      IMPORT c_ptr, c_double
      TYPE(c_ptr), VALUE, INTENT(IN) :: it1, it2
      REAL(c_double), INTENT(OUT) :: t
    END SUBROUTINE util_diff_real_time
  END INTERFACE

END MODULE mo_util_timer
