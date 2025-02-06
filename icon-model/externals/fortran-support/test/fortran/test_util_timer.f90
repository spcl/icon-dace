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

MODULE test_mo_util_timer
  USE FORTUTF

CONTAINS

  SUBROUTINE TEST_cputime
    USE ISO_C_BINDING, ONLY: c_ptr, c_double, c_int
    USE mo_util_timer
    REAL(c_double) :: time1, time2

    CALL TAG_TEST("TEST_gettimeofday")
    time1 = util_gettimeofday()
    time2 = util_gettimeofday()
    CALL ASSERT_ALMOST_EQUAL(time1, time2)
  END SUBROUTINE

END MODULE test_mo_util_timer
