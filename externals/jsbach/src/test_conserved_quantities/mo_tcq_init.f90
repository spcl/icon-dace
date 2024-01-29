!> tcq (test conserved quantities) memory initialisation
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
!>#### Initialization of tcq (test conserved quantities) variables
!>
MODULE mo_tcq_init
#ifndef __NO_JSBACH__

  USE mo_kind,              ONLY: wp

  USE mo_jsb_model_class,   ONLY: t_jsb_model
  USE mo_jsb_tile_class,    ONLY: t_jsb_tile_abstract

  dsl4jsb_Use_processes TCQ_
  dsl4jsb_Use_memory(TCQ_)

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tcq_init

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tcq_init'

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialize tcq process (after memory has been set up)
  !
  SUBROUTINE tcq_init(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':tcq_init'

    dsl4jsb_Def_memory(TCQ_)

    dsl4jsb_Real2D_onDomain :: a_dead_c_tcq
    dsl4jsb_Real2D_onDomain :: a_veg_c_tcq
    dsl4jsb_Real2D_onDomain :: another_dead_c_tcq
    dsl4jsb_Real2D_onDomain :: an_implicit_scaling_tcq
    ! -------------------------------------------------------------------------------------------------- !

    dsl4jsb_Get_memory(TCQ_)

    dsl4jsb_Get_var2D_onDomain(TCQ_, a_dead_c_tcq)
    dsl4jsb_Get_var2D_onDomain(TCQ_, a_veg_c_tcq)
    dsl4jsb_Get_var2D_onDomain(TCQ_, another_dead_c_tcq)
    dsl4jsb_Get_var2D_onDomain(TCQ_, an_implicit_scaling_tcq)

    !>
    !> initialise with a constant depending on lct
    !>
    an_implicit_scaling_tcq(:,:) = 10.0_wp * tile%lcts(1)%lib_id
    a_dead_c_tcq(:,:) = tile%lcts(1)%lib_id
    a_veg_c_tcq(:,:) = 5.0_wp * tile%lcts(1)%lib_id
    another_dead_c_tcq(:,:) = 100.0_wp * tile%lcts(1)%lib_id

  END SUBROUTINE tcq_init

#endif
END MODULE mo_tcq_init
