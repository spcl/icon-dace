!> Contains methods for JSBACH finalization.
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
MODULE mo_jsb_model_final
#ifndef __NO_JSBACH__

  ! USE mo_kind,             ONLY: wp
  ! USE mo_exception,        ONLY: message, finish

  ! USE mo_jsb_class,        ONLY: Get_model
  USE mo_jsb_control,      ONLY: l_timer_host, timer_on !,get_no_of_models, debug_on
  ! USE mo_jsb_model_class,  ONLY: t_jsb_model

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: jsbach_finalize

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_model_final'

CONTAINS

  SUBROUTINE jsbach_finalize

    USE mo_timer, ONLY: print_timer

    ! TYPE(t_jsb_model), POINTER  :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':jsbach_finalize'

    IF (.NOT. l_timer_host .AND. timer_on()) CALL print_timer()

  END SUBROUTINE jsbach_finalize

#endif
END MODULE mo_jsb_model_final
