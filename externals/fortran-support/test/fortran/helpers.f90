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

MODULE helpers
  USE ISO_C_BINDING, ONLY: c_double
  USE mo_io_units, ONLY: find_next_free_unit
  USE mo_exception, ONLY: message

CONTAINS

  ! deterministic sequence of "random" numbers.
  ! (rval in [-1,1], linear congruential generator)
  SUBROUTINE rrand(seed, rval)
    INTEGER, INTENT(INOUT) :: seed
    REAL(c_double), INTENT(OUT) :: rval
    !
    INTEGER, PARAMETER :: rand_m = ISHFT(2, 15) + 1 ! 2**16+1
    INTEGER, PARAMETER :: rand_a = 75
    INTEGER, PARAMETER :: rand_c = 74
    seed = MOD((rand_a*seed + rand_c), rand_m)
    rval = 2.0_c_double*REAL(seed, c_double)/REAL(rand_m - 1, c_double) - 1.0_c_double
  END SUBROUTINE rrand

  SUBROUTINE custom_exit()
    CALL EXIT(1)
  END SUBROUTINE

  SUBROUTINE custom_exit_dummy()
  END SUBROUTINE

  SUBROUTINE open_new_logfile(nerr, file)
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(out) :: nerr

    nerr = find_next_free_unit(10, 20)
    OPEN (unit=nerr, file=file, status='replace', action='write')

  END SUBROUTINE open_new_logfile

  SUBROUTINE open_logfile(nerr, file)
    CHARACTER(len=*), INTENT(in) :: file
    INTEGER, INTENT(in) :: nerr

    OPEN (unit=nerr, file=file, status='old', action='read')

  END SUBROUTINE open_logfile
END MODULE helpers
