!> wlcc (lcc due to wind) initialisation of lcc structure
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
MODULE mo_wlcc_init_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_lcc,             ONLY: init_lcc, count_descendants_for_up_to_2layers, &
    &                               collect_names_of_descendant_leaves_up_to_2layers
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
    &                               BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE

  dsl4jsb_Use_processes WLCC_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: wlcc_init_lcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_wlcc_init_lcc_class'

  INTEGER, PARAMETER :: nr_of_active_cqts = 1
  INTEGER, PARAMETER :: nr_of_passive_cqts = 4

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialise lcc structure for wlcc (lcc due to wind) process
  !
  SUBROUTINE wlcc_init_lcc(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: involved_tiles(:)
    INTEGER :: nr_of_descendants

    INTEGER :: active_cqts(nr_of_active_cqts)
    INTEGER :: passive_cqts(nr_of_passive_cqts)

    CHARACTER(len=*), PARAMETER :: routine = modname//':wlcc_init_lcc'
    ! -------------------------------------------------------------------------------------------------- !

    IF (.NOT. tile%Is_process_calculated(WLCC_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    active_cqts = (/ LIVE_CARBON_CQ_TYPE /)
    passive_cqts = (/ AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE /)

    ! WLCC disturbs area on the leave tiles
    IF (tile%Has_children()) THEN
      CALL count_descendants_for_up_to_2layers(tile, nr_of_descendants)
      ALLOCATE(involved_tiles(nr_of_descendants))
      CALL collect_names_of_descendant_leaves_up_to_2layers(tile, involved_tiles)
    ELSE
      ALLOCATE(involved_tiles(1))
      involved_tiles(1) = tile%name
    END IF

    CALL init_lcc('wlcc', tile, active_cqts, passive_cqts, involved_tiles)

    DEALLOCATE(involved_tiles)

  END SUBROUTINE wlcc_init_lcc

#endif
END MODULE mo_wlcc_init_lcc
