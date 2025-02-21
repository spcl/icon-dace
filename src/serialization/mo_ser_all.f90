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

MODULE mo_ser_all


  IMPLICIT NONE

  PUBLIC :: serialize_all ! This is the only component that has to be available without SERIALIZE

  CONTAINS


  SUBROUTINE serialize_all(nproma, jg, savepoint_base, is_input, opt_lupdate_cpu, opt_id, opt_dt)
    USE mtime, ONLY: datetime

    INTEGER, INTENT(IN) :: nproma, jg
    CHARACTER(LEN=*), INTENT(IN) :: savepoint_base
    LOGICAL, INTENT(IN) :: is_input

    LOGICAL, INTENT(IN), OPTIONAL :: opt_lupdate_cpu
    INTEGER, INTENT(IN), OPTIONAL :: opt_id
    ! use this to pass a datetime that describes the exact current sub-timestep of a nested domain.
    TYPE(datetime), INTENT(IN), OPTIONAL, POINTER :: opt_dt 


  END SUBROUTINE serialize_all

END MODULE mo_ser_all
