!> declare and define indices for bgc material
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
!>#### provides index variables for bgc material compartments and elements in the bgcm_store
!>
MODULE mo_lnd_bgcm_idx

  USE mo_veg_config_class,        ONLY: VEG_PART_LEAF_IDX , VEG_PART_FINE_ROOT_IDX, VEG_PART_COARSE_ROOT_IDX, &
    &                                   VEG_PART_SAP_WOOD_IDX, VEG_PART_HEART_WOOD_IDX, VEG_PART_LABILE_IDX, &
    &                                   VEG_PART_RESERVE_IDX, VEG_PART_FRUIT_IDX, VEG_PART_SEED_BED_IDX
  USE mo_sb_config_class,         ONLY: SB_PART_DOM_IDX , SB_PART_DOM_ASSOC_IDX, SB_PART_SOLUBLE_LITTER_IDX, &
    &                                   SB_PART_POLYMERIC_LITTER_IDX, SB_PART_WOODY_LITTER_IDX, &
    &                                   SB_PART_FUNGI_IDX, SB_PART_MYCORRHIZA_IDX, SB_PART_MICROBIAL_IDX, &
    &                                   SB_PART_RESIDUE_IDX, SB_PART_RESIDUE_ASSOC_IDX
  USE mo_lnd_bgcm_class,          ONLY: ELEM_C_ID, ELEM_N_ID, ELEM_P_ID, ELEM_C13_ID, ELEM_C14_ID, ELEM_N15_ID, &
    &                                   FIRST_ELEM_ID, LAST_ELEM_ID

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ixC, ixN, ixP, ixC13, ixC14, ixN15
  PUBLIC :: FIRST_ELEM_ID, LAST_ELEM_ID
  PUBLIC :: set_bgcm_element_idx
  PUBLIC :: ix_leaf, ix_fine_root, ix_coarse_root, ix_sap_wood, ix_heart_wood, ix_labile, ix_reserve, ix_fruit, ix_seed_bed
  PUBLIC :: ix_dom, ix_dom_assoc, ix_soluable_litter, ix_polymeric_litter, ix_woody_litter, ix_fungi, ix_mycorrhiza, &
    &       ix_microbial, ix_residue, ix_residue_assoc

  ! elements
  INTEGER, SAVE :: ixC   = -1
  INTEGER, SAVE :: ixN   = -1
  INTEGER, SAVE :: ixP   = -1
  INTEGER, SAVE :: ixC13 = -1
  INTEGER, SAVE :: ixC14 = -1
  INTEGER, SAVE :: ixN15 = -1
  ! VEG_
  INTEGER, PARAMETER :: ix_leaf             = VEG_PART_LEAF_IDX
  INTEGER, PARAMETER :: ix_fine_root        = VEG_PART_FINE_ROOT_IDX
  INTEGER, PARAMETER :: ix_coarse_root      = VEG_PART_COARSE_ROOT_IDX
  INTEGER, PARAMETER :: ix_sap_wood         = VEG_PART_SAP_WOOD_IDX
  INTEGER, PARAMETER :: ix_heart_wood       = VEG_PART_HEART_WOOD_IDX
  INTEGER, PARAMETER :: ix_labile           = VEG_PART_LABILE_IDX
  INTEGER, PARAMETER :: ix_reserve          = VEG_PART_RESERVE_IDX
  INTEGER, PARAMETER :: ix_fruit            = VEG_PART_FRUIT_IDX
  INTEGER, PARAMETER :: ix_seed_bed         = VEG_PART_SEED_BED_IDX
  ! SB_
  INTEGER, PARAMETER :: ix_dom              = SB_PART_DOM_IDX
  INTEGER, PARAMETER :: ix_dom_assoc        = SB_PART_DOM_ASSOC_IDX
  INTEGER, PARAMETER :: ix_soluable_litter  = SB_PART_SOLUBLE_LITTER_IDX
  INTEGER, PARAMETER :: ix_polymeric_litter = SB_PART_POLYMERIC_LITTER_IDX
  INTEGER, PARAMETER :: ix_woody_litter     = SB_PART_WOODY_LITTER_IDX
  INTEGER, PARAMETER :: ix_fungi            = SB_PART_FUNGI_IDX
  INTEGER, PARAMETER :: ix_mycorrhiza       = SB_PART_MYCORRHIZA_IDX
  INTEGER, PARAMETER :: ix_microbial        = SB_PART_MICROBIAL_IDX
  INTEGER, PARAMETER :: ix_residue          = SB_PART_RESIDUE_IDX
  INTEGER, PARAMETER :: ix_residue_assoc    = SB_PART_RESIDUE_ASSOC_IDX

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lnd_bgcm_idx'

CONTAINS

  ! ======================================================================================================= !
  !>set the index variables of bgcm element indicess from "model%config%elements_index_map(:)"
  !>
  !>  the "model%config%elements_index_map(:)" is set in mo_jsb_config_class:new_model_config
  !>  based on namelist options
  !>
  SUBROUTINE set_bgcm_element_idx(elements_index_map)
    USE mo_jsb_control,     ONLY: debug_on
    USE mo_exception,       ONLY: message
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER, INTENT(in) :: elements_index_map(:)
    ! ----------------------------------------------------------------------------------------------------- !
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':set_bgcm_element_idx'
    IF (debug_on()) CALL message(TRIM(routine), 'Running ...')

    ixC   = elements_index_map(ELEM_C_ID)
    ixN   = elements_index_map(ELEM_N_ID)
    ixP   = elements_index_map(ELEM_P_ID)
    ixC13 = elements_index_map(ELEM_C13_ID)
    ixC14 = elements_index_map(ELEM_C14_ID)
    ixN15 = elements_index_map(ELEM_N15_ID)
  END SUBROUTINE set_bgcm_element_idx

END MODULE mo_lnd_bgcm_idx
