!> tlcc (test lcc) memory class
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
!>#### Could contain memory definitions for tlcc
!>
!> currently of no use beyond being required by the infrastructure
!>
MODULE mo_tlcc_memory_class
#ifndef __NO_JSBACH__

  USE mo_jsb_model_class,        ONLY: t_jsb_model
  USE mo_jsb_class,              ONLY: Get_model
  USE mo_jsb_memory_class,       ONLY: t_jsb_memory

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_tlcc_memory, max_no_of_vars

  INTEGER, PARAMETER :: max_no_of_vars = 5

  !> Type definition for memory
  TYPE, EXTENDS(t_jsb_memory) :: t_tlcc_memory

  CONTAINS
    PROCEDURE :: Init => Init_tlcc_memory
  END TYPE t_tlcc_memory

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tlcc_memory_class'

CONTAINS

  SUBROUTINE Init_tlcc_memory(mem, prefix, suffix, lct_ids, model_id)

    USE mo_jsb_varlist,       ONLY: t_jsb_varlist
    USE mo_jsb_io,            ONLY: tables
    USE mo_jsb_grid_class,    ONLY: t_jsb_grid, t_jsb_vgrid
    USE mo_jsb_grid,          ONLY: Get_grid, Get_vgrid

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_tlcc_memory), INTENT(inout), TARGET :: mem
    CHARACTER(len=*),    INTENT(in)    :: prefix
    CHARACTER(len=*),    INTENT(in)    :: suffix
    INTEGER,             INTENT(in)    :: lct_ids(:)
    INTEGER,             INTENT(in)    :: model_id
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_model), POINTER :: model
    TYPE(t_jsb_grid),  POINTER :: hgrid                        ! Horizontal grid
    TYPE(t_jsb_vgrid), POINTER :: surface                      ! Vertical grid
    INTEGER :: table

    CHARACTER(len=*), PARAMETER :: routine = modname//':Init_tlcc_memory'
    ! -------------------------------------------------------------------------------------------------- !
    
    model => Get_model(model_id)

    table = tables(1)

    hgrid   => Get_grid(mem%grid_id)
    surface => Get_vgrid('surface')

    !This process currently has no own variables 

  END SUBROUTINE Init_tlcc_memory

#endif
END MODULE mo_tlcc_memory_class
