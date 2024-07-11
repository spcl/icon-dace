!> Initialization of the the atm2land memory.
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
MODULE mo_atm2land_init
#ifndef __NO_QUINCY__

  USE mo_kind,              ONLY: wp
  USE mo_exception,         ONLY: message
  USE mo_jsb_control,       ONLY: debug_on

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: atm2land_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_atm2land_init'

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> Run atm2land init
  !!
  !!
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE atm2land_init(tile)

    USE mo_jsb_class,             ONLY: Get_model
    USE mo_jsb_tile_class,        ONLY: t_jsb_tile_abstract
    USE mo_jsb_model_class,       ONLY: t_jsb_model
    USE mo_jsb_process_class,     ONLY: A2L_
    USE mo_atmland_constants,     ONLY: def_fdiffuse


    ! Use of process configurations (t_PROC_config)
    !dsl4jsb_Use_config(A2L_)

    ! USE PROC_memory
    dsl4jsb_Use_memory(A2L_)


    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    CLASS(t_jsb_tile_abstract), INTENT(inout)     :: tile         !< one tile with data structure for one lct
    ! ---------------------------
    ! 0.2 Local
    TYPE(t_jsb_model),      POINTER              :: model
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':atm2land_init'
    ! ---------------------------
    ! 0.3 Declare Memory
    ! Declare process configuration and memory Pointers
    !dsl4jsb_Def_config(A2L_)
    dsl4jsb_Def_memory(A2L_)
    ! Declare pointers to variables in memory
    !dsl4jsb_Real2D_onDomain      :: fract_par_diffuse
    ! ---------------------------
    ! 0.4 Process Activity, Debug Option
    IF (ASSOCIATED(tile%parent_tile)) RETURN ! should only run on the uppermost tile of the hierarchy

    IF (debug_on()) CALL message(TRIM(routine), 'Setting initial conditions of atm2land memory (quincy) for tile '// &
      &                          TRIM(tile%name))
    ! ---------------------------
    ! 0.5 Get Memory
    model  => Get_model(tile%owner_model_id)
    ! Get process config
    !dsl4jsb_Get_config(A2L_)
    ! Get process memories
    dsl4jsb_Get_memory(A2L_)
    ! Get process variables (Set pointers to variables in memory)
    !dsl4jsb_Get_var2D_onDomain(A2L_, fract_par_diffuse)


    ! ------------------------------------------------------------------------------------------------------------
    ! Go science
    ! ------------------------------------------------------------------------------------------------------------


    !fract_par_diffuse(:,:)               = def_fdiffuse


  END SUBROUTINE atm2land_init

#endif
END MODULE mo_atm2land_init
