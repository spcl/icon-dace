!> helper routines for vegetation (QUINCY)
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
!> For more information on the QUINCY model see: <https://doi.org/10.17871/quincy-model-2019>
!>
!>#### various helper routines for the vegetation process
!>

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_veg_util
#ifndef __NO_JSBACH__

  USE mo_kind,        ONLY: wp
  USE mo_jsb_control, ONLY: debug_on
  USE mo_exception,   ONLY: message

  ! Use of processes in this module
  dsl4jsb_Use_processes VEG_

  ! Use of process memories
  dsl4jsb_Use_memory(VEG_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: reset_veg_fluxes

  CHARACTER(len=*), PARAMETER :: modname = 'mo_veg_util'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> init/reset vegetation fluxes (with/to zero)
  !! 
  !! called prior to all Tasks (pheno_interface)
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE reset_veg_fluxes(tile, options)

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
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':reset_veg_fluxes'

    dsl4jsb_Def_memory(VEG_) 

    dsl4jsb_Real2D_onChunk      :: veg_growth_leaf_carbon
    dsl4jsb_Real2D_onChunk      :: veg_litterfall_leaf_carbon
    dsl4jsb_Real2D_onChunk      :: maint_respiration_pot

    ! Get local variables from options argument
    iblk    = options%iblk
    ics     = options%ics
    ice     = options%ice
    nc      = options%nc

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    dsl4jsb_Get_memory(VEG_) 

    model  => Get_model(tile%owner_model_id)

    dsl4jsb_Get_var2D_onChunk(VEG_, veg_growth_leaf_carbon)           ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, veg_litterfall_leaf_carbon)       ! out
    dsl4jsb_Get_var2D_onChunk(VEG_, maint_respiration_pot)            ! out

    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    ! set zero all flux variables of the assimilation process
    veg_growth_leaf_carbon(:)       = zero
    veg_litterfall_leaf_carbon(:)   = zero
    
    maint_respiration_pot(:)        = zero

    IF (debug_on() .AND. iblk==1) CALL message(TRIM(routine), 'Finished.')

  END SUBROUTINE reset_veg_fluxes

#endif
END MODULE mo_veg_util
