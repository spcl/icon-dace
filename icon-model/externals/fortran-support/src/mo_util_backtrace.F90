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

MODULE mo_util_backtrace

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ftn_util_backtrace

  INTERFACE
    SUBROUTINE util_backtrace() BIND(C)
    END SUBROUTINE util_backtrace
  END INTERFACE

CONTAINS

  SUBROUTINE ftn_util_backtrace()
    CALL util_backtrace()
  END SUBROUTINE ftn_util_backtrace

END MODULE mo_util_backtrace

