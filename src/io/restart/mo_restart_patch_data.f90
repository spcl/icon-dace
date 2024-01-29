! The base CLASS for the PUBLIC INTERFACE for restart writing.
!
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

MODULE mo_restart_patch_data
  USE mo_restart_patch_description, ONLY: t_restart_patch_description
  USE mo_var,                       ONLY: t_var_ptr

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_RestartPatchData

  ! this type stores all the information that we need to know about a patch and its variables
  TYPE, ABSTRACT :: t_RestartPatchData
    TYPE(t_restart_patch_description) :: description
    TYPE(t_var_ptr), ALLOCATABLE :: varData(:)
    INTEGER :: restartType
  CONTAINS
    PROCEDURE(i_construct), DEFERRED :: construct
    PROCEDURE(i_writeData), DEFERRED :: writeData
    PROCEDURE(i_destruct), DEFERRED :: destruct
  END TYPE t_RestartPatchData

  ABSTRACT INTERFACE
  SUBROUTINE i_construct(me, modelType, jg)
    IMPORT t_RestartPatchData
    CLASS(t_RestartPatchData), INTENT(INOUT) :: me
    CHARACTER(*), INTENT(IN) :: modelType
    INTEGER, INTENT(IN) :: jg
  END SUBROUTINE i_construct

  SUBROUTINE i_writeData(me, ncid)
    IMPORT t_RestartPatchData
    CLASS(t_RestartPatchData), INTENT(INOUT), TARGET :: me
    INTEGER, INTENT(IN) :: ncid
  END SUBROUTINE i_writeData

  SUBROUTINE i_destruct(me)
    IMPORT t_RestartPatchData
    CLASS(t_RestartPatchData), INTENT(INOUT) :: me
  END SUBROUTINE i_destruct
  END INTERFACE

END MODULE mo_restart_patch_data
