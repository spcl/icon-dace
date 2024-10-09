!> Summary: provides types and definitions for storing and providing bgc materials
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
!>#### provides types and definitions for storing and providing bgc materials
!>
MODULE mo_lnd_bgcm_store_class
#ifndef __NO_QUINCY__

  USE mo_kind,                ONLY: wp
  USE mo_jsb_pool_class,      ONLY: t_jsb_pool
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p
  USE mo_jsb_varlist,         ONLY: VARNAME_LEN

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lnd_bgcm_store_abstract
  PUBLIC :: SB_BGCM_POOL_ID, SB_BGCM_DELTA_POOLS_ID, SB_BGCM_FORMATION_ID, SB_BGCM_LOSS_ID, &
    &       SB_BGCM_TRANSPORT_ID, SB_BGCM_MYCO_EXPORT_ID
  PUBLIC :: VEG_BGCM_POOL_ID, VEG_BGCM_TOTAL_BIO_ID, VEG_BGCM_GROWTH_ID,               &
    &       VEG_BGCM_LITTERFALL_ID, VEG_BGCM_FRAC_ALLOC_ID,                            &
    &       VEG_BGCM_EXUDATION_ID, VEG_BGCM_ESTABLISHMENT_ID, VEG_BGCM_RESERVE_USE_ID, &
    &       VEG_BGCM_PP_FUEL_ID, VEG_BGCM_PP_PAPER_ID, VEG_BGCM_PP_FIBERBOARD_ID,      &
    &       VEG_BGCM_PP_OIRW_ID, VEG_BGCM_PP_PV_ID, VEG_BGCM_PP_SAWNWOOD_ID, VEG_BGCM_FPROD_DECAY_ID

  ! ======================================================================================================= !
  !>
  !> IDs of all bgc materials that are to be used in the store
  ENUM, BIND(C)
    ENUMERATOR ::                   &
      & SB_BGCM_POOL_ID = 1,        &   !< soil biogeochemistry bgcm pool
      & SB_BGCM_DELTA_POOLS_ID,     &   !< soil biogeochemistry bgcm delta pools
      & SB_BGCM_FORMATION_ID,       &   !< soil biogeochemistry bgcm formation
      & SB_BGCM_LOSS_ID,            &   !< soil biogeochemistry bgcm loss
      & SB_BGCM_TRANSPORT_ID,       &   !< soil biogeochemistry bgcm transport
      & SB_BGCM_MYCO_EXPORT_ID,     &   !< soil biogeochemistry bgcm mycorrhiza export
      & VEG_BGCM_POOL_ID,           &   !< vegetation bgcm pool
      & VEG_BGCM_TOTAL_BIO_ID,      &   !< vegetation bgcm total biomass
      & VEG_BGCM_GROWTH_ID,         &   !< vegetation bgcm growth
      & VEG_BGCM_LITTERFALL_ID,     &   !< vegetation bgcm litterfall
      & VEG_BGCM_FRAC_ALLOC_ID,     &   !< vegetation bgcm frac_alloc
      & VEG_BGCM_EXUDATION_ID,      &   !< vegetation bgcm exudation
      & VEG_BGCM_ESTABLISHMENT_ID,  &   !< vegetation bgcm establishment
      & VEG_BGCM_RESERVE_USE_ID,    &   !< vegetation bgcm reserve_use
      & VEG_BGCM_PP_FUEL_ID,        &   !< vegetation bgcm fuel product pool
      & VEG_BGCM_PP_PAPER_ID,       &   !< vegetation bgcm paper product pool
      & VEG_BGCM_PP_FIBERBOARD_ID,  &   !< vegetation bgcm fiberboard product pool
      & VEG_BGCM_PP_OIRW_ID,        &   !< vegetation bgcm other industrial roundwood product pool
      & VEG_BGCM_PP_PV_ID,          &   !< vegetation bgcm plywood and veneer product pool
      & VEG_BGCM_PP_SAWNWOOD_ID,    &   !< vegetation bgcm sawnwood product pool
      & VEG_BGCM_FPROD_DECAY_ID         !< vegetation bgcm flux from product pool decay
  END ENUM

  ! ======================================================================================================= !
  !>
  !> Type used to collect pointer for a certain 1l bgc material (i.e. only elements, no compartements)
  !>
  TYPE :: t_lnd_1l_bgc_material_collection
    INTEGER                        :: nr_of_elems   !< number of elements (array size)
    TYPE(t_jsb_var_p), ALLOCATABLE :: collection_1l_p(:) !< pointer array for access to the vars
    REAL(wp),              POINTER :: mt_1l_2d_bgcm(:,:,:)   !< if 2d var (elements, nc, nblks)
    REAL(wp),              POINTER :: mt_1l_3d_bgcm(:,:,:,:) !< if 3d var (elements, nc, nsoil, nblks)
  END TYPE t_lnd_1l_bgc_material_collection

  ! ======================================================================================================= !
  !>
  !> Type used to collect pointer for a certain 2l bgc material (i.e. with compartments)
  !>
  TYPE :: t_lnd_2l_bgc_material_collection
    INTEGER                        :: nr_of_parts   !< number of compartments (size dim 1)
    INTEGER                        :: nr_of_elems   !< number of elements (size dim 2)
    TYPE(t_jsb_var_p), ALLOCATABLE :: collection_2l_p(:,:)   !< pointer array for access to the vars
    REAL(wp),              POINTER :: mt_2l_2d_bgcm(:,:,:,:) !< if 2d var (compartments, elements, nc, nblks)
    REAL(wp),              POINTER :: mt_2l_3d_bgcm(:,:,:,:,:)
      !< if 3d var (compartments, elements, nc, nsoil, nblks)
  END TYPE t_lnd_2l_bgc_material_collection

  ! ======================================================================================================= !
  !>
  !> Type used on the tile for bgc material bookkeeping and data storage on matrices
  !>
  TYPE, ABSTRACT :: t_lnd_bgcm_store_abstract
    CHARACTER(len=VARNAME_LEN), ALLOCATABLE :: bgc_material_name(:) !< name (len: nr of bgc materials)
    INTEGER, ALLOCATABLE :: bgc_material_ID(:)  !< ID - ENUM (len: same ++ also below)
    INTEGER, ALLOCATABLE :: nr_of_levels(:)     !< number of levels (currently: 1 or 2)
    INTEGER, ALLOCATABLE :: nr_of_dimensions(:) !< number of dimensions (currently: 2 or 3)
    INTEGER, ALLOCATABLE :: idx_in_store(:)     !< index in the according store
    INTEGER :: nr_of_bgc_materials       = 0    !< total number of bgc materials in the store
    INTEGER :: nr_of_1l_2d_bgc_materials = 0    !< integer to keep track of 1l 2d bgc materials
    INTEGER :: nr_of_1l_3d_bgc_materials = 0    !< integer to keep track of 1l 3d bgc materials
    INTEGER :: nr_of_2l_2d_bgc_materials = 0    !< integer to keep track of 2l 2d bgc materials
    INTEGER :: nr_of_2l_3d_bgc_materials = 0    !< integer to keep track of 2l 3d bgc materials
    TYPE(t_lnd_1l_bgc_material_collection), ALLOCATABLE :: store_1l_2d_bgcms(:)
      !< bcg material collection for one layer bgc materials with elements with two dimensions
    TYPE(t_lnd_1l_bgc_material_collection), ALLOCATABLE :: store_1l_3d_bgcms(:)
      !< bcg material collection for one layer bgc materials with elements with three dimensions
    TYPE(t_lnd_2l_bgc_material_collection), ALLOCATABLE :: store_2l_2d_bgcms(:)
      !< bcg material collection for bgc materials with compartments with elements with two dimensions
    TYPE(t_lnd_2l_bgc_material_collection), ALLOCATABLE :: store_2l_3d_bgcms(:)
      !< bcg material collection for bgc materials with compartments with elements with three dimensions
    CONTAINS
      PROCEDURE(Init_interface),         DEFERRED :: Init
      PROCEDURE(Append_bgcms_interface), DEFERRED :: Append_bgcms
  END TYPE t_lnd_bgcm_store_abstract

  ABSTRACT INTERFACE
    SUBROUTINE Init_interface(this)
      IMPORT :: t_lnd_bgcm_store_abstract
      CLASS(t_lnd_bgcm_store_abstract), INTENT(inout) :: this
    END SUBROUTINE Init_interface
  END INTERFACE

  ABSTRACT INTERFACE
    SUBROUTINE Append_bgcms_interface(this, main_bgcm, nblks, nproma, nsoil, tile_name, idx_bgcm, &
        &                             idx_1l_2d_bgcm, idx_2l_2d_bgcm, idx_1l_3d_bgcm, idx_2l_3d_bgcm)
      IMPORT :: t_lnd_bgcm_store_abstract, t_jsb_pool
      CLASS(t_lnd_bgcm_store_abstract), INTENT(inout) :: this
      CLASS(t_jsb_pool),                INTENT(in)    :: main_bgcm
      INTEGER,                          INTENT(in)    :: nblks, nproma, nsoil
      CHARACTER(len=*),                 INTENT(in)    :: tile_name
      INTEGER,                          INTENT(inout) :: idx_bgcm
      INTEGER,                          INTENT(inout) :: idx_1l_2d_bgcm, idx_2l_2d_bgcm
      INTEGER,                          INTENT(inout) :: idx_1l_3d_bgcm, idx_2l_3d_bgcm

    END SUBROUTINE Append_bgcms_interface
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lnd_bgcm_store_class'

CONTAINS

#endif
END MODULE mo_lnd_bgcm_store_class
