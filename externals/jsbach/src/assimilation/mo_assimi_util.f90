!> diverse subroutines and functions for the assimilation process (QUINCY)
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

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_assimi_util
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

  ! Use of processes in this module
  dsl4jsb_Use_processes ASSIMI_

  ! Use of process memories
  dsl4jsb_Use_memory(ASSIMI_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reset_assimi_fluxes    ! may be called each timestep prior to all Tasks

  CHARACTER(len=*), PARAMETER :: modname = 'mo_assimi_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> init/reset assimilation fluxes (with/to zero)
  !! 
  !! called prior to all Tasks (assimi_interface)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_assimi_fluxes(tile, options)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_task_class,        ONLY: t_jsb_task_options
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_math_constants,    ONLY: zero
    
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    TYPE(t_jsb_task_options),   INTENT(in)        :: options      !< model options

    ! Locals
    TYPE(t_jsb_model),      POINTER       :: model
    INTEGER                               :: iblk, ics, ice, nc
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_assimi_fluxes'

    dsl4jsb_Def_memory(ASSIMI_) 

    ! ASSIMI_ 2D
    dsl4jsb_Real2D_onChunk      :: gross_assimilation
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C13
    dsl4jsb_Real2D_onChunk      :: gross_assimilation_C14
    dsl4jsb_Real2D_onChunk      :: net_assimilation
    dsl4jsb_Real2D_onChunk      :: net_assimilation_boc
    dsl4jsb_Real2D_onChunk      :: maint_respiration_leaf
    ! ASSIMI_ 3D
    dsl4jsb_Real3D_onChunk      :: gross_assimilation_cl
    dsl4jsb_Real3D_onChunk      :: net_assimilation_cl
    dsl4jsb_Real3D_onChunk      :: maint_respiration_leaf_cl

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (.NOT. tile%Is_process_active(ASSIMI_)) RETURN

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    ! Get process memories
    dsl4jsb_Get_memory(ASSIMI_) 

    ! Get process variables (Set pointers to variables in memory)
    model  => Get_model(tile%owner_model_id)
    ! ---------------------------
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_C13)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  gross_assimilation_C14)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  net_assimilation)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  net_assimilation_boc)
    dsl4jsb_Get_var2D_onChunk(ASSIMI_,  maint_respiration_leaf)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,  gross_assimilation_cl)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,  net_assimilation_cl)
    dsl4jsb_Get_var3D_onChunk(ASSIMI_,  maint_respiration_leaf_cl)



    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! set zero all flux variables of the assimilation process
    gross_assimilation(:)           = zero
    gross_assimilation_C13(:)       = zero
    gross_assimilation_C14(:)       = zero
    net_assimilation(:)             = zero
    net_assimilation_boc(:)         = zero
    maint_respiration_leaf(:)       = zero

    gross_assimilation_cl(:,:)      = zero
    net_assimilation_cl(:,:)        = zero
    maint_respiration_leaf_cl(:,:)  = zero

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE reset_assimi_fluxes

#endif
END MODULE mo_assimi_util
