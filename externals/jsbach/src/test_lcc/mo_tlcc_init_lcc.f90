!> tlcc (test lcc) lcc structure initialisation
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
!>#### Initialization of the tlcc lcc structure
!>
MODULE mo_tlcc_init_lcc
#ifndef __NO_JSBACH__

  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_lcc,             ONLY: init_lcc
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, DEAD_CARBON_CQ_TYPE

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: tlcc_init_lcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_tlcc_init_lcc_class'
  CHARACTER(len=*), PARAMETER :: procname = 'tlcc'
  
  INTEGER, PARAMETER :: nr_of_active_cqts = 1 
  INTEGER, PARAMETER :: nr_of_passive_cqts = 1
  INTEGER, PARAMETER :: nr_of_involved_descendant_tiles = 5

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialise lcc structure for tlcc (test lcc) process
  !
  SUBROUTINE tlcc_init_lcc(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: active_cqts(nr_of_active_cqts)
    INTEGER :: passive_cqts(nr_of_passive_cqts)
    CHARACTER(len=SHORT_NAME_LEN) :: involved_child_tiles(nr_of_involved_descendant_tiles)
    ! -------------------------------------------------------------------------------------------------- !
    
    ! JN-TODO: reconsider what to define here and how
    IF (.NOT. tile%name .EQ. 'veg') RETURN

    active_cqts = (/ LIVE_CARBON_CQ_TYPE /)
    passive_cqts = (/ DEAD_CARBON_CQ_TYPE /)

    involved_child_tiles = [character(len=SHORT_NAME_LEN) :: 'pft03', 'pft04' , 'pft05', 'pft06', 'pft11' ]

    CALL init_lcc(procname, tile, active_cqts, passive_cqts, involved_child_tiles)

  END SUBROUTINE tlcc_init_lcc

#endif
END MODULE mo_tlcc_init_lcc
