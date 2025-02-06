!> Contains basic model structure
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
MODULE mo_jsb_class
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class, ONLY: t_jsb_model, t_jsb_model_m

  IMPLICIT NONE
  PRIVATE

  ! t_jsb needs to be public for the NEC compiler
  PUBLIC :: t_jsb, jsbach, get_model

  TYPE t_jsb
     INTEGER :: no_of_models
     TYPE(t_jsb_model_m), POINTER :: models(:)
  END TYPE t_jsb

  TYPE(t_jsb) :: jsbach

  INTERFACE get_model
     MODULE PROCEDURE get_model_by_id
     !MODULE PROCEDURE get_model_by_name
  END INTERFACE get_model

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_class'

CONTAINS

!!$  FUNCTION get_grid_by_id(id) RESULT(grid)
!!$
!!$    USE mo_jsb_grid_class, ONLY: t_jsb_grid
!!$
!!$    INTEGER, INTENT(in)        :: id
!!$    TYPE(t_jsb_grid), POINTER :: grid
!!$
!!$    grid => jsbach%models(id)%m%grid
!!$
!!$  END FUNCTION get_grid_by_id

  FUNCTION get_model_by_id(id) RESULT(model)

    INTEGER, INTENT(in)        :: id
    TYPE(t_jsb_model), POINTER :: model

    CHARACTER(len=*), PARAMETER :: routine = modname//':get_model_by_id'

    model => jsbach%models(id)%m

  END FUNCTION get_model_by_id

#endif
END MODULE mo_jsb_class
