!> fage (forest age) initialisation of lcc structure
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
!>#### Initialization of the lcc structure for fage (forest age)
!>
MODULE mo_fage_init_lcc
#ifndef __NO_JSBACH__

  USE mo_exception,           ONLY: message
  USE mo_jsb_control,         ONLY: debug_on
  USE mo_jsb_tile_class,      ONLY: t_jsb_tile_abstract
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN
  USE mo_jsb_lcc,             ONLY: init_lcc
  USE mo_jsb_cqt_class,       ONLY: LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, &
    &                               BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE

  dsl4jsb_Use_processes FAGE_

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: fage_init_lcc

  CHARACTER(len=*), PARAMETER :: modname = 'mo_fage_init_lcc_class'
  CHARACTER(len=*), PARAMETER :: procname = 'fage'

  INTEGER, PARAMETER :: nr_of_active_cqts = 0   !< fage currently does not have any active cqts
  INTEGER, PARAMETER :: nr_of_passive_cqts = 5

CONTAINS

  ! ====================================================================================================== !
  !
  !> Initialise lcc structure for fage (forest age) process
  !
  SUBROUTINE fage_init_lcc(tile)

    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
    ! -------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = modname//':fage_init_lcc'

    INTEGER :: active_cqts(nr_of_active_cqts)
    INTEGER :: passive_cqts(nr_of_passive_cqts)
    INTEGER :: nr_of_involved_tiles

    CHARACTER(len=SHORT_NAME_LEN), ALLOCATABLE :: involved_child_tiles(:)
    ! -------------------------------------------------------------------------------------------------- !

    IF (.NOT. tile%Is_process_calculated(FAGE_)) RETURN

    IF (debug_on()) CALL message(TRIM(routine), 'Starting on tile '//TRIM(tile%name)//' ...')

    passive_cqts = (/ LIVE_CARBON_CQ_TYPE, AG_DEAD_C_CQ_TYPE, BG_DEAD_C_CQ_TYPE, PRODUCT_CARBON_CQ_TYPE, FLUX_C_CQ_TYPE /)

    nr_of_involved_tiles = tile%Get_no_of_children()
    ALLOCATE(involved_child_tiles(nr_of_involved_tiles))
    CALL tile%Get_children_names(involved_child_tiles)

    CALL init_lcc(procname, tile, active_cqts, passive_cqts, involved_child_tiles)

    DEALLOCATE(involved_child_tiles)

  END SUBROUTINE fage_init_lcc

#endif
END MODULE mo_fage_init_lcc
