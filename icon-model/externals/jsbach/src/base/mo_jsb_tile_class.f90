!> Contains abstract tile class that contains the surface structure and memory state
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
!>### Definition of abstract class and methods for tiles
!> The most important components of a tile instance are
!>
!> - the fraction of a grid box that the tile occupies
!> - the memories of all processes (i.e. variables) owned by this tile
!> - a list of LCT (land cover type) id's that are represented on this tile
!> - a store pointing to the element variables of all BGC materials across all processes
!> - a list of actions to be applied on this tile for the different processes
!> - a list of aggregators used to aggregate variables of the child tiles to this tile

!NEC$ options "-finline-file=externals/jsbach/src/base/mo_jsb_control.pp-jsb.f90"

MODULE mo_jsb_tile_class
#ifndef __NO_JSBACH__

  USE mo_jsb_control,         ONLY: debug_on
  USE mo_kind,                ONLY: wp, dp
  USE mo_exception,           ONLY: message, finish
  USE mo_io_units,            ONLY: filename_max
  USE mo_jsb_impl_constants,  ONLY: SHORT_NAME_LEN

  USE mo_hsm_class,           ONLY: t_Hsm, t_State, t_Message
  USE mo_jsb_memory_class,    ONLY: t_jsb_memory_p, t_jsb_memory
  USE mo_jsb_lct_class,       ONLY: t_jsb_lct, &
    &                               LAND_TYPE, VEG_TYPE, BARE_TYPE, GLACIER_TYPE, LAKE_TYPE
  USE mo_jsb_var_class,       ONLY: t_jsb_var_p, t_jsb_var, t_jsb_var_real2d, t_jsb_var_real3d
  USE mo_jsb_cqt_class,       ONLY: t_jsb_consQuan_p
  USE mo_jsb_lcc_class,       ONLY: t_jsb_lcc_proc_p
#ifndef __NO_QUINCY__
  USE mo_lnd_bgcm_store,      ONLY: t_lnd_bgcm_store
#endif

#ifdef _OPENACC
  USE openacc
#define __acc_attach(ptr) CALL acc_attach(ptr)
#else
#define __acc_attach(ptr)
#endif

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_jsb_tile_abstract
  PUBLIC :: t_jsb_aggregator, t_jsb_aggregator_p, t_jsb_aggregator_weighted_by_fract

  ! Declarations for aggregators

  !> Abstract type for aggregator operators
  TYPE, ABSTRACT :: t_jsb_aggregator
    CHARACTER(len=25) :: name    !< Aggregators are distinguished by using unique names
  CONTAINS
    PROCEDURE(aggregate_interface_2d), DEFERRED             :: Aggregate_2d !< Aggregate operator for 2d variables
    PROCEDURE(aggregate_interface_3d), DEFERRED             :: Aggregate_3d !< Aggregate operator for 3d variables
    GENERIC                                                 :: Aggregate => Aggregate_2d, Aggregate_3d
  END TYPE t_jsb_aggregator

  !> A type to hold vectors of pointers to t_jsb_aggregator objects
  TYPE t_jsb_aggregator_p
    CLASS(t_jsb_aggregator), POINTER :: p => NULL()
  END TYPE t_jsb_aggregator_p

  ! Declarations for tiles
  !
  TYPE :: real_2d_p
    REAL(wp), POINTER :: p(:,:) => NULL()
  END TYPE
  TYPE :: real_3d_p
    REAL(wp), POINTER :: p(:,:,:) => NULL()
  END TYPE
  !
  TYPE, ABSTRACT, EXTENDS(t_State) :: t_jsb_tile_abstract
    TYPE(t_jsb_lct),        ALLOCATABLE    :: lcts(:) !< List of land cover types for tile;
                                                      !< the first entry is always the primary LCT for this tile;
                                                      !< the remaining lcts are those inherited for aggregation
    ! Convinience settings to easily check primary LCT types and those of the child tiles
    LOGICAL                                :: is_vegetation = .FALSE.       !< whether primary lct is VEG_TYPE
    LOGICAL                                :: is_bare = .FALSE.             !< whether primary lct is BARE_TYPE
    LOGICAL                                :: is_land = .FALSE.             !< whether primary lct is LAND_TYPE
    LOGICAL                                :: is_glacier = .FALSE.          !< whether primary lct is GLACIER_TYPE
    LOGICAL                                :: is_lake = .FALSE.             !< whether primary lct is LAKE_TYPE
    LOGICAL                                :: contains_vegetation = .FALSE. !< whether any lct of tile is VEG_TYPE
    LOGICAL                                :: contains_bare = .FALSE.       !< whether any lct of tile is BARE_TYPE
    LOGICAL                                :: contains_land = .FALSE.       !< whether any lct of tile is LAND_TYPE
    LOGICAL                                :: contains_glacier = .FALSE.    !< whether any lct of tile is GLACIER_TYPE
    LOGICAL                                :: contains_lake = .FALSE.       !< whether any lct of tile is LAKE_TYPE
    LOGICAL                                :: contains_soil = .FALSE.       !< whether any lct of tile is one that has soil
    !
    INTEGER,                  ALLOCATABLE :: process_action(:)        !< List of `action` for all processes
    REAL(wp),                     POINTER :: fract     (:,:)          !< Fraction of this tile (relative to grid box)
    REAL(wp),                     POINTER :: fract_old (:,:)          !< Old fractions of this tile (relative to grid box)
    CHARACTER(LEN=filename_max)           :: fract_filename           !< Name of file to read tile fractions from. The variable
                                                                      !< names must be either "fract[_max]_"//tile%name,
                                                                      !< or specified by fract_varname (in the use case).
    CHARACTER(LEN=20)                     :: fract_varname                 !< Name of variable containing fractions for this tile
    LOGICAL,                  ALLOCATABLE :: l_fract_children_changed(:,:) !< indicates whether the fraction is changed
    TYPE(t_jsb_aggregator_p), ALLOCATABLE :: aggregators(:)                !< List of possible aggregators for this tile
    !
    INTEGER                               :: nr_of_cqts = 0                !< number of conserved quantity types in use on the tile
    TYPE(t_jsb_consQuan_p),   ALLOCATABLE :: conserved_quantities(:)       !< list of conserved quantities sort by type
    TYPE(t_jsb_lcc_proc_p),   ALLOCATABLE :: lcc_processes(:)              !< list of lcc structures indexed by proc id
#ifndef __NO_QUINCY__
    TYPE(t_lnd_bgcm_store),       POINTER :: bgcm_store          => NULL() !< the bgcm store collecting pools and their elements on this tile
#endif
    TYPE(t_jsb_memory_p),         POINTER :: mem(:)                        !< Memory for all processes on the tile
    !
    CLASS(t_jsb_tile_abstract),   POINTER :: self_tile           => NULL() !< pointer to the tile to enable assigning it as target
    ! The next three pointers basically repeat what's already in the base type, but they save "SELECT TYPE"s
    CLASS(t_jsb_tile_abstract),   POINTER :: parent_tile         => NULL() !< pointer to the parent file
    CLASS(t_jsb_tile_abstract),   POINTER :: first_child_tile    => NULL() !< pointer to the first child tile
    CLASS(t_jsb_tile_abstract),   POINTER :: next_sibling_tile   => NULL() !< pointer to the next sibling tile
    !
    TYPE(real_2d_p),              POINTER :: ptrs2d_cache(:,:,:) => NULL() !< cache of pointers to 2d vars used for GPU porting
    TYPE(real_3d_p),              POINTER :: ptrs3d_cache(:,:,:) => NULL() !< cache of pointers to 3d vars used for GPU porting
    !
    INTEGER                               :: grid_id = -1                  !< ID of the grid used for the tile
    INTEGER                               :: owner_model_id = 0            !< ID of model instance (needed e.g. for nesting)

  CONTAINS
    PROCEDURE :: Add_lct
    PROCEDURE :: Add_lct_to_parents
    PROCEDURE :: Get_parent_memory
    PROCEDURE :: Set_process_action => Set_process_action_tile
    PROCEDURE :: Set_process_action_parents
    PROCEDURE :: Register_aggregator
    PROCEDURE :: Get_aggregator
    PROCEDURE :: Set_fraction_1d                                  !< Set fractions on tile (chunk)
    PROCEDURE :: Set_fraction_2d                                  !< Set fractions on tile (domain)
    GENERIC   :: Set_fraction => Set_fraction_1d, Set_fraction_2d !< Set fractions on tile
    PROCEDURE :: Get_fraction_1d
    PROCEDURE :: Get_fraction_2d
    GENERIC   :: Get_fraction => Get_fraction_1d, Get_fraction_2d
    PROCEDURE :: Set_fract_ptr
    PROCEDURE :: Set_fract_old_ptr
    PROCEDURE :: Get_first_child_tile
    PROCEDURE :: Get_next_sibling_tile
    PROCEDURE :: Get_prev_sibling_tile
    PROCEDURE :: Get_children_names
    PROCEDURE :: Has_conserved_quantities
    PROCEDURE :: Cache_GPU_pointers
    !
    PROCEDURE(Init_interface),           DEFERRED :: Init
    PROCEDURE(Init_vars_interface),      DEFERRED :: Init_vars
#ifndef __NO_QUINCY__
    PROCEDURE(Init_vars_interface),      DEFERRED :: Count_and_classify_bgc_materials
    PROCEDURE(Init_vars_interface),      DEFERRED :: Collect_bgc_materials
#endif
    PROCEDURE(Init_vars_interface),      DEFERRED :: Count_conserved_quantities
    PROCEDURE(Init_vars_interface),      DEFERRED :: Collect_conserved_variables
    PROCEDURE(Init_fractions_interface), DEFERRED :: Init_fractions
    PROCEDURE(Print_interface),          DEFERRED :: Print
    PROCEDURE(Handler_interface),        DEFERRED :: Handler
    PROCEDURE(Check_on_process),         DEFERRED :: Is_process_active
    PROCEDURE(Check_on_process),         DEFERRED :: Is_process_calculated
    PROCEDURE(Check_on_process),         DEFERRED :: Has_process_memory
    PROCEDURE                                     :: Is_last_process_tile
  END type t_jsb_tile_abstract

  ! Deferred subroutine interfaces (must be provided as part of t_jsb_tile_abstract).
  ! The subroutines themselves are located in the jsb_tile module.
  ABSTRACT INTERFACE
    !> Interface for deferred subroutine `Init`
    SUBROUTINE Init_interface(this, varlist_name, prefix, suffix, grid_id, in_var_groups)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
      CHARACTER(len=*),           INTENT(in)    :: varlist_name
      CHARACTER(len=*),           INTENT(in)    :: prefix
      CHARACTER(len=*),           INTENT(in)    :: suffix
      INTEGER,                    INTENT(in)    :: grid_id
      LOGICAL,                    INTENT(in)    :: in_var_groups
    END SUBROUTINE Init_interface

    !> Interface for deferred subroutine `Init_fractions`
    SUBROUTINE Init_fractions_interface(this, varlist_name, prefix, suffix, l_fixed_fractions, l_rel_fractions)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
      CHARACTER(len=*),           INTENT(in)    :: varlist_name
      CHARACTER(len=*),           INTENT(in)    :: prefix
      CHARACTER(len=*),           INTENT(in)    :: suffix
      LOGICAL,                    INTENT(in)    :: l_fixed_fractions
      LOGICAL,                    INTENT(in)    :: l_rel_fractions
    END SUBROUTINE Init_fractions_interface

    !> Interface for varous deferred subroutines
    SUBROUTINE Init_vars_interface(this)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    END SUBROUTINE Init_vars_interface

    !> Interface for various deferred functions
    LOGICAL FUNCTION Check_on_process(this, iproc)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(in) :: this
      INTEGER, INTENT(in) :: iproc
    END FUNCTION Check_on_process

    !> Interface for deferred function `Handler`
    FUNCTION Handler_interface(this, msg_in) RESULT(return_ptr)
      IMPORT :: t_jsb_tile_abstract, t_Hsm, t_Message
      CLASS(t_jsb_tile_abstract),  INTENT(inout) :: this
      CLASS(t_Message),            INTENT(in)    :: msg_in
      CLASS(t_Message),            POINTER       :: return_ptr
    END FUNCTION Handler_interface

    !> Interface for deferred subroutine `Print`
    SUBROUTINE Print_interface(this)
      IMPORT :: t_jsb_tile_abstract
      CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    END SUBROUTINE Print_interface

  END INTERFACE

  ! Interfaces for deferred subroutines of abstract t_jsb_aggregator type

  ABSTRACT INTERFACE
    !> Interface for deferred subroutine `Aggregate_2d` of `t_jsb_aggregator`
    SUBROUTINE aggregate_interface_2d(this, tile, var, ics, ice, iblk, caller)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real2d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this            !< Aggregator
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real2d),     INTENT(inout) :: var             !< Variable with 2 dimensions
      INTEGER,                    INTENT(in)    :: ics, ice, iblk  !< Array indices of current chunk
      CHARACTER(len=*),           INTENT(in)    :: caller          !< For debugging: Where was aggregator called?
    END SUBROUTINE aggregate_interface_2d

    !> Interface for deferred subroutine `Aggregate_3d` of `t_jsb_aggregator`
    SUBROUTINE aggregate_interface_3d(this, tile, var, ics, ice, iblk,caller)
      IMPORT :: t_jsb_aggregator, t_jsb_tile_abstract, t_jsb_var_real3d
      CLASS(t_jsb_aggregator),    INTENT(inout) :: this
      CLASS(t_jsb_tile_abstract), INTENT(inout) :: tile
      TYPE(t_jsb_var_real3d),     INTENT(inout) :: var             !< Variable with 3 dimensions
      INTEGER,                    INTENT(in)    :: ics, ice, iblk  !< Array indices of current chunk
      CHARACTER(len=*),           INTENT(in)    :: caller          !< For debugging: Where was aggregator called?
    END SUBROUTINE aggregate_interface_3d

  END INTERFACE

  !> Implementation of a concrete aggregator for aggregation by weighted fractions
  TYPE, EXTENDS(t_jsb_aggregator) :: t_jsb_aggregator_weighted_by_fract
    !> Variables needed to calculate aggregation weights
    REAL(wp), ALLOCATABLE :: fractions(:,:,:) !< Weighting fraction for each child of tile @domain (nproma, nblock, no_children)
    REAL(wp), ALLOCATABLE :: fract_sum(:,:)   !< Sum of fractions for fraction weighted averaging
    REAL(wp), ALLOCATABLE :: data_sum(:,:)    !< Sum of data*fraction for fraction weighted averaging
  CONTAINS
    !> Functions to make use of the aggregator
    PROCEDURE :: Set_fractions_1d => Set_fractions_aggregator_1d  !< Set fraction field for chunk. TODO: Rename "1d" to "chunk"
    PROCEDURE :: Set_fractions_2d => Set_fractions_aggregator_2d  !< Set fraction field for entire domain. TODO: "2d"=>"domain"
    GENERIC   :: Set_fractions    => Set_fractions_1d, Set_fractions_2d !< Generic function to set tile fractions
    PROCEDURE :: Aggregate_2d     => Aggregate_weighted_by_fract_2d     !< Specific implementation of agg. of 2d vars
    PROCEDURE :: Aggregate_3d     => Aggregate_weighted_by_fract_3d     !< Specific implementation of agg. of 3d vars
    FINAL     :: Finalize_aggregator_weighted_by_fract            !< Destructor to deallocate fields when aggregator is deallocated
  END TYPE t_jsb_aggregator_weighted_by_fract

  INTERFACE t_jsb_aggregator_weighted_by_fract
    PROCEDURE Create_aggregator_weighted_by_fract                 !< Constructor of aggregator (call to initialize)
  END INTERFACE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_jsb_tile_class'

CONTAINS

  FUNCTION Get_parent_memory(this, iproc) RESULT(return_ptr)

    CLASS(t_jsb_tile_abstract),  INTENT(in) :: this
    INTEGER,                     INTENT(in) :: iproc
    CLASS(t_jsb_memory),         POINTER   :: return_ptr

    CLASS(*), POINTER :: parent

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_parent_memory'

    parent => this%Get_parent()
    IF (ASSOCIATED(parent)) THEN
      SELECT TYPE (parent)
      CLASS IS (t_jsb_tile_abstract)
        return_ptr => parent%mem(iproc)%p
      END SELECT
    ELSE
      return_ptr => NULL()
    END IF

  END FUNCTION Get_parent_memory

  SUBROUTINE Set_process_action_tile(this, iproc, action, action_to_parents)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    INTEGER,                    INTENT(in)    :: iproc
    INTEGER,                    INTENT(in)    :: action
    INTEGER, OPTIONAL,          INTENT(in)    :: action_to_parents

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_process_action_tile'

    this%process_action(iproc) = action

    ! Propagate action action_to_parents to all parents of this tile
    IF (PRESENT(action_to_parents)) CALL this%Set_process_action_parents(iproc, action_to_parents)

  END SUBROUTINE Set_process_action_tile

  SUBROUTINE Set_process_action_parents(this, iproc, action)

    CLASS(t_jsb_tile_abstract), INTENT(inout), TARGET :: this
    INTEGER,                    INTENT(in)            :: iproc
    INTEGER,                    INTENT(in)            :: action

    CLASS(t_jsb_tile_abstract), POINTER :: tile

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_process_action_parents'

    tile => this
    DO WHILE (ASSOCIATED(tile%parent))
      SELECT TYPE (parent=>tile%parent)
      CLASS IS (t_jsb_tile_abstract)
        parent%process_action(iproc) = action
        tile => parent
      END SELECT
    END DO

  END SUBROUTINE Set_process_action_parents

  ! ====================================================================================================== !
  !
  !> adds given lct to this tile and sets the "is" and/or "contains" properties on the tile accordingly
  !
  SUBROUTINE Add_lct(this, lct)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this  !< this tile
    TYPE(t_jsb_lct),            INTENT(in)    :: lct   !< the landcover type to add for this tile
    ! -------------------------------------------------------------------------------------------------- !
    TYPE(t_jsb_lct), ALLOCATABLE :: temp_lcts(:)
    INTEGER :: ilct, nlct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_lct'
    ! -------------------------------------------------------------------------------------------------- !
    IF (.NOT. ALLOCATED(this%lcts)) THEN

      ALLOCATE(this%lcts(1))
      this%lcts(1) = lct

      ! the "is" properties are only set for the primary lct
      IF (ASSOCIATED(this%parent)) THEN    ! Root tile has no primary LCT
        SELECT CASE (lct%id)
        CASE (VEG_TYPE)
          this%is_vegetation = .TRUE.
        CASE (BARE_TYPE)
          this%is_bare = .TRUE.
        CASE (LAND_TYPE)
          this%is_land = .TRUE.
        CASE (GLACIER_TYPE)
          this%is_glacier = .TRUE.
        CASE (LAKE_TYPE)
          this%is_lake = .TRUE.
        END SELECT
      END IF

    ELSE

      nlct = SIZE(this%lcts)
      DO ilct=1,nlct
        ! If this lct is already registered in the tile%lcts, do nothing and return
        IF (this%lcts(ilct)%id == lct%id) RETURN
      END DO
      ALLOCATE(temp_lcts(nlct+1))
      temp_lcts(1:nlct) = this%lcts
      temp_lcts(nlct+1) = lct
      CALL MOVE_ALLOC(temp_lcts, this%lcts)

    END IF !(.NOT. ALLOCATED(this%lcts))

    ! the "contains" properties are set according to all lcts of this tile
    SELECT CASE (lct%id)
    CASE (VEG_TYPE)
      this%contains_vegetation = .TRUE.
    CASE (BARE_TYPE)
      this%contains_bare = .TRUE.
    CASE (LAND_TYPE)
      this%contains_land = .TRUE.
    CASE (GLACIER_TYPE)
      this%contains_glacier = .TRUE.
    CASE (LAKE_TYPE)
      this%contains_lake = .TRUE.
    END SELECT

    this%contains_soil       = this%contains_vegetation .OR. this%contains_bare .OR. this%contains_land

  END SUBROUTINE Add_lct

  ! ====================================================================================================== !
  !
  !> adds given lct to all parent tiles of this tile
  !
  SUBROUTINE Add_lct_to_parents(this, lct)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this  !< tile whose parent is to be considered
    TYPE(t_jsb_lct),            INTENT(in)    :: lct   !< the landcover type to add to the parent tiles
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_State), POINTER :: current

    CHARACTER(len=*), PARAMETER :: routine = modname//':Add_lct_to_parents'
    ! -------------------------------------------------------------------------------------------------- !
    current => this%parent

    DO WHILE (ASSOCIATED(current))
      SELECT TYPE (current)
      CLASS IS (t_jsb_tile_abstract)
        CALL current%Add_lct(lct)
      END SELECT
      current => current%parent
    END DO

  END SUBROUTINE Add_lct_to_parents

  SUBROUTINE Set_fraction_1d(this, ics, ice, iblk, fract, fract_old)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    INTEGER,                    INTENT(in)    :: ics, ice, iblk
    REAL(wp), OPTIONAL,         INTENT(in)    :: &
      & fract(:),        &
      & fract_old(:)

    INTEGER :: ilct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fraction_1d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract) /= ice-ics+1) &
       & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_old(ics:ice, iblk) = this%fract(ics:ice, iblk)
      this%fract    (ics:ice, iblk) = fract(:)
      !$ACC UPDATE DEVICE(this%fract(ics:ice,iblk), this%fract_old(ics:ice,iblk)) ASYNC(1)

      IF (ASSOCIATED(this%parent)) THEN
        SELECT TYPE (parent=>this%parent)
        CLASS IS (t_jsb_tile_abstract)
          WHERE (this%fract(ics:ice, iblk) /= this%fract_old(ics:ice, iblk))
            parent%l_fract_children_changed(ics:ice, iblk) = .TRUE.
          END WHERE
          IF (ANY(parent%l_fract_children_changed(ics:ice,iblk))) THEN
            DO ilct=1,SIZE(parent%lcts)
              IF (this%lcts(1)%id == parent%lcts(ilct)%id) THEN
                WHERE (parent%l_fract_children_changed(ics:ice,iblk))
                  parent%lcts(ilct)%fract(ics:ice,iblk) = &
                    & parent%lcts(ilct)%fract(ics:ice,iblk) - this%fract_old(ics:ice,iblk) + this%fract(ics:ice,iblk)
                END WHERE
                EXIT
              END IF
            END DO
          END IF
        END SELECT
      END IF
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old) /= ice-ics+1) &
       & CALL finish(TRIM(routine), 'Dimension mismatch')

      this%fract_old(ics:ice, iblk) = fract_old(:)
      !$ACC UPDATE DEVICE(this%fract_old(ics:ice,iblk)) ASYNC(1)
    END IF

  END SUBROUTINE Set_fraction_1d

  SUBROUTINE Set_fraction_2d(this, fract, fract_old)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), OPTIONAL,         INTENT(in)    :: &
      & fract(:,:),        &
      & fract_old(:,:)

    INTEGER :: ilct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fraction_2d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract,1) /= SIZE(this%fract,1) .OR. SIZE(fract,2) /= SIZE(this%fract,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      IF (ANY(fract(:,:) /= this%fract(:,:))) THEN

        this%fract_old(:,:) = this%fract(:,:)
        this%fract    (:,:) = fract(:,:)
        !$ACC UPDATE DEVICE(this%fract, this%fract_old) ASYNC(1)

        IF (ASSOCIATED(this%parent)) THEN
          SELECT TYPE (parent=>this%parent)
          CLASS IS (t_jsb_tile_abstract)
            WHERE (this%fract(:,:) /= this%fract_old(:,:))
              parent%l_fract_children_changed(:,:) = .TRUE.
            END WHERE
            IF (ANY(parent%l_fract_children_changed(:,:))) THEN
              DO ilct=1,SIZE(parent%lcts)
                IF (this%lcts(1)%id == parent%lcts(ilct)%id) THEN
                  WHERE (parent%l_fract_children_changed(:,:))
                    parent%lcts(ilct)%fract(:,:) = parent%lcts(ilct)%fract(:,:) - this%fract_old(:,:) + this%fract(:,:)
                  END WHERE
                  EXIT
                END IF
              END DO
            END IF
          END SELECT
        END IF

      END IF
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old,1) /= SIZE(this%fract_old,1) .OR. SIZE(fract_old,2) /= SIZE(this%fract_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')

      IF (ANY(fract_old(:,:) /= this%fract_old(:,:))) THEN
        this%fract_old(:,:) = fract_old(:,:)
        !$ACC UPDATE DEVICE(this%fract_old) ASYNC(1)
      END IF

    END IF

  END SUBROUTINE Set_fraction_2d

  SUBROUTINE Get_fraction_1d(this, ics, ice, iblk, fract, fract_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    INTEGER,                    INTENT(in)    :: ics, ice, iblk
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:),        &
      & fract_old(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_fraction_1d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract(:) = this%fract(ics:ice, iblk)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old) /= ice-ics+1) CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:) = this%fract_old(ics:ice, iblk)
    END IF

  END SUBROUTINE Get_fraction_1d

  SUBROUTINE Get_fraction_2d(this, fract, fract_old)

    CLASS(t_jsb_tile_abstract), INTENT(in)    :: this
    REAL(wp), OPTIONAL,         INTENT(inout) :: &
      & fract(:,:),        &
      & fract_old(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_fraction_2d'

    IF (PRESENT(fract)) THEN
      IF (SIZE(fract,1) /= SIZE(this%fract,1) .OR. SIZE(fract,2) /= SIZE(this%fract,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract(:,:) = this%fract(:,:)
    END IF

    IF (PRESENT(fract_old)) THEN
      IF (SIZE(fract_old,1) /= SIZE(this%fract_old,1) .OR. SIZE(fract_old,2) /= SIZE(this%fract_old,2)) &
        & CALL finish(TRIM(routine), 'Dimension mismatch')
      fract_old(:,:) = this%fract_old(:,:)
    END IF

  END SUBROUTINE Get_fraction_2d

  ! ======================================================================================================= !
  !>
  !> Sets the pointer on the tile to the passed pointer AND also updates the GPU pointer address
  !>
  ! This is an example why the encapsulation within the setter is important: not only the fraction
  ! on the tile needs to be set but also the GPU pointer
  SUBROUTINE Set_fract_ptr(this, ptr)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), POINTER                         :: ptr(:,:)
    ! -------------------------------------------------------------------------------------------------- !
    this%fract => ptr
    __acc_attach(this%fract)

  END SUBROUTINE Set_fract_ptr

  SUBROUTINE Set_fract_old_ptr(this, ptr)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this
    REAL(wp), POINTER                         :: ptr(:,:)

    this%fract_old => ptr
    __acc_attach(this%fract_old)

  END SUBROUTINE Set_fract_old_ptr

  ! ======================================================================================================= !
  !>
  !> Registers a new aggregator for use on this tile
  !>
  SUBROUTINE Register_aggregator(this, aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(inout)        :: this       !< Tile to which to associate new aggregator
    CLASS(t_jsb_aggregator),    INTENT(in),   TARGET :: aggregator !< Aggregator to associate

    TYPE(t_jsb_aggregator_p), ALLOCATABLE :: temp_aggregators(:)
    INTEGER :: n
    CHARACTER(len=*), PARAMETER :: routine = modname//':Register_aggregator'

    ! Extend tile's list of aggregators by one entry
    IF (ALLOCATED(this%aggregators)) THEN
      n = SIZE(this%aggregators)
    ELSE
      n = 0
    END IF
    ALLOCATE(temp_aggregators(n+1))
    IF (ALLOCATED(this%aggregators)) temp_aggregators(1:n) = this%aggregators ! Copy all old aggregators to temporaray array
    temp_aggregators(n+1) = t_jsb_aggregator_p(aggregator)                    ! Put new aggregator in the new entry of the list
    CALL move_ALLOC(temp_aggregators, this%aggregators)               ! Replace list with temporary and deallocate temporary

    IF (debug_on()) CALL message(TRIM(routine), 'Registered aggregator '//TRIM(aggregator%name)//' for tile '//TRIM(this%name))

  END SUBROUTINE Register_aggregator

  ! ======================================================================================================= !
  !>
  !> Returns an aggregator specified by it's name in the tile's aggregator list
  !> If one suitable is found either a pointer to it (default) or a copy of it is returned
  !>
  FUNCTION Get_aggregator(this, name, copy) RESULT(aggregator)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this        !< Tile to obtain valid aggregator for
    CHARACTER(len=*),           INTENT(in) :: name        !< Name of requested aggregator
    LOGICAL, OPTIONAL,          INTENT(in) :: copy        !< Return a copy instead of a pointer (default: .FALSE.)
    CLASS(t_jsb_aggregator),    POINTER    :: aggregator  !< Returned aggregator (NULL if no aggregator with given name found)

    INTEGER :: i
    LOGICAL :: l_copy

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_aggregator'

    l_copy = .FALSE.
    IF (PRESENT(copy)) l_copy = copy

    aggregator => NULL()
    DO i=1,SIZE(this%aggregators)
      IF (TRIM(this%aggregators(i)%p%name) == TRIM(name)) THEN
        IF (l_copy) THEN
          ALLOCATE(aggregator, source=this%aggregators(i)%p)
        ELSE
          aggregator => this%aggregators(i)%p
        END IF
        EXIT
      END IF
    END DO

    IF (.NOT.(ASSOCIATED(aggregator))) &
      & CALL finish(TRIM(routine), 'Aggregator "'//TRIM(name)//'" not found on tile '//TRIM(this%name))

  END FUNCTION Get_aggregator

  ! ======================================================================================================= !
  !>
  !> Constructs the "weighted_by_fract" aggregator by allocating its internal fields
  !>
  !> Uses the size of the domain and the number of children of the tile
  !>
  FUNCTION Create_aggregator_weighted_by_fract(grid_id, no_of_children) RESULT(aggregator)

    USE mo_jsb_grid_class, ONLY: t_jsb_grid
    USE mo_jsb_grid,       ONLY: Get_grid

    INTEGER, INTENT(in) :: grid_id        !< Id of horizontal grid to get domain size from
    INTEGER, INTENT(in) :: no_of_children !< number of child tiles which are to be aggregated on this tile
    TYPE(t_jsb_aggregator_weighted_by_fract), POINTER    :: aggregator     !< Constructed "weighted_by_fract" aggregator

    TYPE(t_jsb_grid), POINTER :: grid

    CHARACTER(len=*), PARAMETER :: routine = modname//':Create_aggregator_weighted_by_fract'

    IF (no_of_children < 1) CALL finish(TRIM(routine), 'This should not happen!')

    ALLOCATE(aggregator)
    aggregator%name = 'weighted_by_fract'

    grid => Get_grid(grid_id)
    ALLOCATE(aggregator%fractions(grid%nproma, grid%nblks, no_of_children))
    ALLOCATE(aggregator%data_sum (grid%nproma, grid%nblks))
    ALLOCATE(aggregator%fract_sum(grid%nproma, grid%nblks))
    !$ACC ENTER DATA CREATE(aggregator, aggregator%fractions, aggregator%data_sum, aggregator%fract_sum)
  END FUNCTION Create_aggregator_weighted_by_fract

  ! ======================================================================================================= !
  !>
  !> Destructs an "weighted_by_fract" aggregator by deallocating its internal fields
  !>
  SUBROUTINE Finalize_aggregator_weighted_by_fract(aggregator)

    TYPE(t_jsb_aggregator_weighted_by_fract) :: aggregator !< Aggregator to destruct

    CHARACTER(len=*), PARAMETER :: routine = modname//':Finalize_aggregator_weighted_by_fract'

    !CALL message(TRIM(routine), 'Finalizing aggregator '//TRIM(aggregator%name))

    IF (ALLOCATED(aggregator%fractions)) THEN
      !$ACC EXIT DATA DELETE(aggregator%fractions)
      DEALLOCATE(aggregator%fractions)
    END IF
    IF (ALLOCATED(aggregator%data_sum)) THEN
      !$ACC EXIT DATA DELETE(aggregator%data_sum)
      DEALLOCATE(aggregator%data_sum)
    END IF
    IF (ALLOCATED(aggregator%fract_sum)) THEN
      !$ACC EXIT DATA DELETE(aggregator%fract_sum)
      DEALLOCATE(aggregator%fract_sum)
    END IF

  END SUBROUTINE Finalize_aggregator_weighted_by_fract

  ! ====================================================================================================== !
  !
  !> Updates the part of the fraction field of an aggregator belonging to the specified chunk for all child tiles
  !
  ! TODO: Change name to Set_fractions_aggregator_chunk
  SUBROUTINE Set_fractions_aggregator_1d(this, tile, ics, ice, iblk)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this           !< Aggregator to set fractions for
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile           !< "Parent tile" of aggregator
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk !< Array indices of chunk

    INTEGER :: no_children, i
    CLASS(t_jsb_tile_abstract), POINTER :: current
    ! REAL(wp), ALLOCATABLE :: fract(:), fract_sum(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fractions_aggregator_1d'

    ! Check if child tile fractions have changed and, if necessary, update fractions of aggregator
    ! by looping over children of tile
    IF (ANY(tile%l_fract_children_changed(ics:ice,iblk))) THEN
      no_children = tile%Get_no_of_children()
      current => tile%first_child_tile
      DO i=1,no_children
        CALL current%Get_fraction(ics, ice, iblk, fract=this%fractions(ics:ice,iblk,i))
        current => current%next_sibling_tile
      END DO

      !$ACC UPDATE DEVICE(this%fractions(ics:ice,iblk,:)) ASYNC(1)

      !! Section for debugging
      !! It needs to be adapted to absolute tile fractions!
      !!   => replace 1._wp by parent%fract
      !
      ! ALLOCATE(fract(ice-ics+1))
      ! ALLOCATE(fract_sum(ice-ics+1))
      ! fract(:) = tile%fract(ics:ice,iblk)
      ! fract_sum(:) = SUM(this%fractions(ics:ice,iblk,:), DIM=2)
      ! ! IF (ASSOCIATED(tile%parent) .AND. ANY( fract_sum(:) < 1._wp - EPSILON(1._wp) .AND. fract(:) > 0._wp)) THEN
      ! !   CALL message(TRIM(routine), 'Sum of child tile fractions on tile '//TRIM(tile%name)//': '//&
      ! !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
      ! !     &          all_print=.TRUE.)
      ! !   CALL finish(TRIM(routine), 'Sum of tile fractions is not 1 on tile')
      ! ! END IF
      ! DEALLOCATE(fract, fract_sum)

      tile%l_fract_children_changed(ics:ice,iblk) = .FALSE.
    END IF

  END SUBROUTINE Set_fractions_aggregator_1d

  ! ====================================================================================================== !
  !
  !> Updates the fraction field of an aggregator for the entire domain for all child tiles
  !
  ! TODO: Change name to Set_fractions_aggregator_domain
  SUBROUTINE Set_fractions_aggregator_2d(this, tile)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this !< Aggregator to set fractions for
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile !< Tile of this aggregator instance

    INTEGER :: no_children, i
    CLASS(t_jsb_tile_abstract), POINTER :: current
    ! REAL(wp), ALLOCATABLE :: fract(:), fract_sum(:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Set_fractions_aggregator_2d'

    ! print*,'EEE ',iblk,ics,ice,ANY(tile%l_fract_children_changed(ics:ice,iblk))
    ! Check if child tile fractions have changed and, if necessary, update fractions of aggregator
    ! by looping over children of tile
    IF (ANY(tile%l_fract_children_changed(:,:))) THEN
      no_children = tile%Get_no_of_children()
      current => tile%first_child_tile
      DO i=1,no_children
        CALL current%Get_fraction_2d(fract=this%fractions(:,:,i))
        current => current%next_sibling_tile
      END DO

      !$ACC UPDATE DEVICE(this%fractions(:,:,:)) ASYNC(1)

      ! ALLOCATE(fract(ice-ics+1))
      ! ALLOCATE(fract_sum(ice-ics+1))
      ! fract(:) = tile%fract(ics:ice,iblk)
      ! fract_sum(:) = SUM(this%fractions(ics:ice,iblk,:), DIM=2)
      ! ! IF (ASSOCIATED(tile%parent) .AND. ANY( fract_sum(:) < 1._wp - EPSILON(1._wp) .AND. fract(:) > 0._wp)) THEN
      ! !   CALL message(TRIM(routine), 'Sum of child tile fractions on tile '//TRIM(tile%name)//': '//&
      ! !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
      ! !     &          all_print=.TRUE.)
      ! !   CALL finish(TRIM(routine), 'Sum of tile fractions is not 1 on tile')
      ! ! END IF
      ! DEALLOCATE(fract, fract_sum)

      tile%l_fract_children_changed(:,:) = .FALSE.
    END IF

  END SUBROUTINE Set_fractions_aggregator_2d

  ! ====================================================================================================== !
  !>
  !> Aggregation of a 2d variable from the child tiles by weighted fractions
  !>
  SUBROUTINE Aggregate_weighted_by_fract_2d(this, tile, var, ics, ice, iblk, caller)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this           !< Aggregator
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile           !< Tile of this aggregator instance
    TYPE(t_jsb_var_real2d),                    INTENT(inout) :: var            !< to be aggregated 2d variable
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk !< chunk info
    CHARACTER(len=*),                          INTENT(in)    :: caller         !< calling aggregation routine

    REAL(wp), POINTER :: ptr(:,:)
    TYPE(real_2d_p), POINTER :: tile_ptrs2d_cache(:,:,:)

    INTEGER :: no_children, icount,ichild, i, iproc, j
    LOGICAL :: l_aggregate_all

    CHARACTER(len=256) :: routine

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    IF (debug_on()) THEN
      routine = modname//':Aggregate_weighted_by_fract_2d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = modname//':Aggregate_weighted_by_fract_2d'
    END IF

    no_children = tile%Get_no_of_children()
    IF (no_children == 0) CALL finish(routine, 'This should never happen!')

    ! Set the fractions of the child tiles (on the current chunk)
    CALL this%Set_fractions(tile, ics, ice, iblk)

    ! Assure the current variable belongs to a process
    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    IF (.NOT. ASSOCIATED(tile%mem(iproc)%p)) RETURN

    ! Check, that the bookkeeping for this variable on its child tiles has been conducted,
    ! and that number of children matches the child indices.
    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    ! Needed for GPUs
    tile_ptrs2d_cache => tile%ptrs2d_cache
    l_aggregate_all = var%l_aggregate_all

    ! print*, routine//": working on var ", var%full_name

    icount = 0  ! Counter: How many children have been considered
    ichild = 0  ! Index of the child
    IF (no_children > 1) THEN
      DO i=1,no_children
        IF ( var%child_idx(i) > 0 ) THEN
          icount = icount + 1
          ichild = i
        ELSE IF (l_aggregate_all) THEN
          icount = icount + 1
          ichild = i
        END IF
      END DO
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(this, tile_ptrs2d_cache, var) ASYNC(1)

    IF (no_children == 1 .AND. var%child_idx(1) > 0) THEN
      ! If there is only one child nothing needs to be aggregated ...
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j = ics, ice
        ! ... and the variable on the current tile equals the var on the child tile.
        var%ptr(j,iblk) = tile_ptrs2d_cache(1,var%child_idx(1),iproc)%p(j,iblk)
      END DO
    ELSE IF (no_children > 1) THEN
      ! With more then one child tile we need to do the weighted averaging.
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j = ics, ice
        this%data_sum (j,iblk) = 0._wp
        this%fract_sum(j,iblk) = 0._wp
      END DO

      ! Calculate the (weighted) sum of the data and the sum up the fractions
      !
      !$ACC LOOP SEQ
      DO i=1,no_children
        ! Check if the var exists on this tile (could be a tile on which the according process is not running)
        IF ( var%child_idx(i) > 0 ) THEN
          ! The variable exists on the child tile
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO j = ics, ice
            ! Multilply variable of the child with the weight and add it to the data sum
            this%data_sum(j,iblk) = this%data_sum(j,iblk) &
              & + this%fractions(j,iblk,i) * tile_ptrs2d_cache(i,var%child_idx(i),iproc)%p(j,iblk)
            ! -> The cache is a huge array of pointers with i iterating over the number of children
            !    (second dim: number of variables, third dim: number of processes).
            ! => Without gpu we would use
            !      tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p%ptr2d(j,iblk)

            ! Add weight to the weight sum
            this%fract_sum(j,iblk) = this%fract_sum(j,iblk) + this%fractions(j,iblk,i)
          END DO
        ELSE IF (l_aggregate_all) THEN
          ! The variable is not defined on the this child tile, but we still want to consider its fraction
          ! (e.g. for fluxes in m-2, that are not calculated on all child tiles; same as adding value zero).
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO j = ics, ice
            this%fract_sum(j,iblk) = this%fract_sum(j,iblk) + this%fractions(j,iblk,i)
          END DO
        END IF
      END DO

      ! print*,'AAA ',ics,ice,MINVAL(this%data_sum(ics:ice)), MAXVAL(this%data_sum(ics:ice)),&
        ! & MINVAL(this%fract_sum(ics:ice)), MAXVAL(this%fract_sum(ics:ice))

      ! IF (ANY(fract_sum(:) > 1._wp + EPSILON(1._wp))) THEN
      !   ! print*, fract_sum
      !   CALL message(TRIM(routine), 'On tile '//TRIM(tile%name)//': fract_sum > 1 - '// &
      !     &                         real2string(MINVAL(fract_sum(:)))//' - '//real2string(MAXVAL(fract_sum(:))), &
      !     &          all_print=.TRUE.)
      !   CALL finish(TRIM(routine),  '    ... this should not happen')
      ! END IF

      ! Calculate the weighted average
      !
      IF (icount == 1) THEN  ! Only one child tile needs to be considered.
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO j = ics, ice
          IF (this%fract_sum(j,iblk) > 0._wp) THEN
            ! ichild: the only child that contained the variable
            var%ptr(j,iblk) = tile_ptrs2d_cache(ichild,var%child_idx(ichild),iproc)%p(j,iblk)
          END IF
        END DO
      ELSE IF (icount > 1) THEN  ! Several child tiles need to be considered.
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO j = ics, ice
          IF (this%fract_sum(j,iblk) > 0._wp) THEN
            ! Divide the data sum by the sum of the fractions
            var%ptr(j,iblk) = this%data_sum(j,iblk) / this%fract_sum(j,iblk)
          ! ELSE
          !   IF (TRIM(tile%name) == 'box') var%ptr(j,iblk) = -1.e4
          END IF
        END DO
      ELSE
        STOP 'Aggregate_weighted_by_fract_2d: This should not happen (icount=0)!'
      END IF

    END IF

    !$ACC END PARALLEL

    IF (l_aggregate_all) THEN
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//'; all_tiles)')
    ELSE
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')
    END IF

  END SUBROUTINE Aggregate_weighted_by_fract_2d


  ! ====================================================================================================== !
  !>
  !> Aggregation of a 3d variable from the child tiles by weighted fractions
  !>
  SUBROUTINE Aggregate_weighted_by_fract_3d(this, tile, var, ics, ice, iblk, caller)

    CLASS(t_jsb_aggregator_weighted_by_fract), INTENT(inout) :: this           !< Aggregator
    CLASS(t_jsb_tile_abstract),                INTENT(inout) :: tile           !< Tile of this aggregator instance
    TYPE(t_jsb_var_real3d),                    INTENT(inout) :: var            !< to be aggregated 3d variable
    INTEGER,                                   INTENT(in)    :: ics, ice, iblk !< chunk info
    CHARACTER(len=*),                          INTENT(in)    :: caller         !< calling aggregation routine

    CLASS(t_jsb_var), POINTER :: ptr
    TYPE(real_3d_p), POINTER :: tile_ptrs3d_cache(:,:,:)

    INTEGER :: no_children, icount, ichild, i, iproc, ilev, nlev, j
    LOGICAL :: l_aggregate_all

    CHARACTER(len=256) :: routine

    IF (.NOT. ASSOCIATED(var%ptr)) RETURN

    IF (debug_on()) THEN
      routine = modname//':Aggregate_weighted_by_fract_3d'//'(called by '//TRIM(caller)//')'
    ELSE
      routine = modname//':Aggregate_weighted_by_fract_3d'
    END IF

    no_children = tile%Get_no_of_children()
    IF (no_children == 0) CALL finish(routine, 'This should never happen!')

    ! Set the fractions of the child tiles (on the current chunk)
    CALL this%Set_fractions(tile, ics, ice, iblk)

    ! Assure the current variable belongs to a process
    iproc = var%owner_proc_id
    IF (iproc < 0) CALL finish(TRIM(routine), 'Unknown process for variable '//TRIM(var%name)//' on tile '//TRIM(tile%name))

    IF (.NOT. ASSOCIATED(tile%mem(iproc)%p)) RETURN

    ! Check, that the bookkeeping for this variable on its child tiles has been conducted,
    ! and that number of children matches the child indices.
    IF (.NOT. ALLOCATED(var%child_idx)) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' not allocated')
    ELSE IF (SIZE(var%child_idx) /= no_children) THEN
      CALL finish(TRIM(routine), 'ERROR - child_idx of '//TRIM(var%name)//' on tile '//TRIM(tile%name)//' has wrong size')
    END IF

    ! Needed for GPUs
    tile_ptrs3d_cache => tile%ptrs3d_cache ! the 3d var pointer cache on this tile
    l_aggregate_all = var%l_aggregate_all

    nlev = SIZE(var%ptr,2)

    ! print*, routine//": working on var ", var%full_name

    icount = 0
    ichild = 0
    IF (no_children > 1) THEN
      DO i=1,no_children
        IF ( var%child_idx(i) > 0 ) THEN
          icount = icount + 1
          ichild = i
        ELSE IF (l_aggregate_all) THEN
          icount = icount + 1
          ichild = i
        END IF
      END DO
    END IF

    !$ACC PARALLEL DEFAULT(PRESENT) PRESENT(this, tile_ptrs3d_cache, var) ASYNC(1)

    IF (no_children == 1 .AND. var%child_idx(1) > 0) THEN
      ! If there is only one child nothing needs to be aggregated ...

      !$ACC LOOP SEQ
      DO ilev = 1, nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO j = ics, ice
          ! ... and the variable on the current tile equals the var on the child tile.
          var%ptr(j,ilev,iblk) = tile_ptrs3d_cache(1,var%child_idx(1),iproc)%p(j,ilev,iblk)
        END DO
      END DO

    ELSE IF (no_children > 1) THEN
      ! With more then one child tile we need to do the weighted averaging.
      !$ACC LOOP GANG(STATIC: 1) VECTOR
      DO j = ics, ice
        this%fract_sum(j,iblk) = 0._wp
      END DO
      !$ACC LOOP SEQ
      DO i=1,no_children
        ! Check if the var exists on this tile (could be a tile on which the according process is not running)
        IF (var%child_idx(i) > 0) THEN
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO j = ics, ice
            this%fract_sum(j,iblk) = this%fract_sum(j,iblk) + this%fractions(j,iblk,i)
          END DO
        ELSE IF (l_aggregate_all) THEN
          ! The variable is not defined on the this child tile, but we still want to consider its fraction
          ! (e.g. for fluxes in m-2, that are not calculated on all child tiles; same as adding value zero).
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO j = ics, ice
            ! Add weight for this child i to the weight sum
            this%fract_sum(j,iblk) = this%fract_sum(j,iblk) + this%fractions(j,iblk,i)
          END DO
        END IF
      END DO

      !$ACC LOOP SEQ
      DO ilev = 1, nlev
        !$ACC LOOP GANG(STATIC: 1) VECTOR
        DO j = ics, ice
          this%data_sum (j,iblk) = 0._wp
        END DO

        !$ACC LOOP SEQ
        DO i=1,no_children
          IF (var%child_idx(i) > 0) THEN
            !$ACC LOOP GANG(STATIC: 1) VECTOR
            DO j = ics, ice
              ! Multilply variable of the child with the weight and add it to the data sum
              this%data_sum(j,iblk) = this%data_sum(j,iblk) &
                & + this%fractions(j,iblk,i) * tile_ptrs3d_cache(i,var%child_idx(i),iproc)%p(j,ilev,iblk)
              ! -> The cache is a huge array of pointers with i iterating over the number of children.
              !    second dim: number of variables, third dim: number of processes
              ! => Without gpu we would use
              !      tile%mem(iproc)%p%children(i)%p%vars(var%child_idx(i))%p%ptr3d(j,ilev,iblk)
            END DO
          END IF
        END DO

        ! Calculate the weighted average
        !
        IF (icount == 1) THEN ! Only one child tile needs to be considered
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO j = ics, ice
            IF (this%fract_sum(j,iblk) > 0._wp) THEN
              ! ichild: the only child that contained the variable
              var%ptr(j,ilev,iblk) = tile_ptrs3d_cache(ichild,var%child_idx(ichild),iproc)%p(j,ilev,iblk)
            END IF
          END DO
        ELSE IF (icount > 1) THEN ! Several child tiles need to be considered
          !$ACC LOOP GANG(STATIC: 1) VECTOR
          DO j = ics, ice
            IF (this%fract_sum(j,iblk) > 0._wp) THEN
              ! Divide the data sum by the sum of the fractions
              var%ptr(j,ilev,iblk) = this%data_sum(j,iblk) / this%fract_sum(j,iblk)
            END IF
          END DO
        ELSE
          STOP 'Aggregate_weighted_by_fract_3d: This should not happen (icount=0)!'
        END IF
      END DO

    END IF

    !$ACC END PARALLEL

    IF (l_aggregate_all) THEN
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//'; all_tiles)')
    ELSE
      IF (debug_on() .AND. iblk==1) CALL message('     ... ', TRIM(var%name)//' ('//TRIM(this%name)//')')
    END IF

  END SUBROUTINE Aggregate_weighted_by_fract_3d

  FUNCTION Get_first_child_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_tile_abstract), INTENT(in)  :: this
    CLASS(t_jsb_tile_abstract), POINTER     :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_first_child_tile'

    next => this%first_child
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find first child tile')

  END FUNCTION Get_first_child_tile

  FUNCTION Get_next_sibling_tile(this) RESULT(return_ptr)

    CLASS(t_jsb_tile_abstract), INTENT(in)  :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr

    CLASS(t_State), POINTER :: next

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_tile_class:Get_next_sibling_tile'

    next => this%next_sibling
    IF (.NOT. ASSOCIATED(next)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (next)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => next
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(next)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find next sibling tile')

  END FUNCTION Get_next_sibling_tile

  ! ====================================================================================================== !
  !
  !> Returns the previous sibling of this tile, if it has a previous sibling
  !
  FUNCTION Get_prev_sibling_tile(this) RESULT(return_ptr)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(in)  :: this
    CLASS(t_jsb_tile_abstract), POINTER :: return_ptr
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_State), POINTER :: prev

    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_tile_class:Get_prev_sibling_tile'
    ! -------------------------------------------------------------------------------------------------- !

    prev => this%prev_sibling
    IF (.NOT. ASSOCIATED(prev)) THEN
      return_ptr => NULL()
      RETURN
    END IF

    SELECT TYPE (prev)
    CLASS IS (t_jsb_tile_abstract)
      return_ptr => prev
    CLASS IS (t_State)
      CALL finish(TRIM(routine), 'current tile is of type t_State, should be t_jsb_tile_abstract')
    CLASS DEFAULT
      CALL finish(TRIM(routine), 'Unkown type for tile')
    END SELECT

    NULLIFY(prev)

    IF (.NOT. ASSOCIATED(return_ptr)) &
      & CALL finish(TRIM(routine), 'Could not find previous sibling tile')

  END FUNCTION Get_prev_sibling_tile

  ! ====================================================================================================== !
  !
  !> Collect the names of the children of this tile
  !
  SUBROUTINE Get_children_names(this, children_names)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_tile_abstract), INTENT(inout)  :: this
    CHARACTER(len=SHORT_NAME_LEN), INTENT(inout) :: children_names(:)
    ! -------------------------------------------------------------------------------------------------- !
    INTEGER :: i_tile, nr_of_children
    CLASS(t_jsb_tile_abstract),  POINTER :: current_tile
    CHARACTER(len=*), PARAMETER :: routine = 'mo_jsb_tile_class:Get_children_names'
    ! -------------------------------------------------------------------------------------------------- !

    nr_of_children = this%Get_no_of_children()

    IF (nr_of_children /= size(children_names)) &
      & CALL finish(routine, 'Call with invalid array size (size does not equal nr of children).')

    current_tile => this%Get_first_child_tile()
    i_tile = 0
    DO WHILE (ASSOCIATED(current_tile))
      i_tile = i_tile + 1
      children_names(i_tile) = current_tile%name
      current_tile => current_tile%Get_next_sibling_tile()
    ENDDO

  END SUBROUTINE Get_children_names


  ! If conserved quantities were collected for this tile
  FUNCTION Has_conserved_quantities(this) RESULT(return_value)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    LOGICAL                   :: return_value

    return_value = this%nr_of_cqts /= 0

  END FUNCTION Has_conserved_quantities

  LOGICAL FUNCTION Is_last_process_tile(this, process_id)

    CLASS(t_jsb_tile_abstract), INTENT(in) :: this
    INTEGER,                    INTENT(in) :: process_id

    Is_last_process_tile = .FALSE.

    IF (this%Has_children()) THEN
      ! TODO INHERIT_ == 1 ; can't use INHERIT_ from mo_jsb_process_class because of cyclic dependency
      IF (this%first_child_tile%process_action(process_id) == 1 .AND. .NOT. ASSOCIATED(this%next_sibling)) THEN
        Is_last_process_tile = .TRUE.
      END IF
    ELSE
      Is_last_process_tile = this%Is_last_leaf()
    END IF

  END FUNCTION Is_last_process_tile

  ! Cache pointers to variables for aggregation (in order to get it to work on GPUs).
  ! Note: This is called during initialization from model_init for each tile that is not a leaf. If var
  !       pointers get reallocated from time to time, the second part of this routine must be called from
  !       within the time loop. But be mindful of OpenMP threads ... do this only once on each processor!
  !
  SUBROUTINE Cache_GPU_pointers(this)

    CLASS(t_jsb_tile_abstract), INTENT(inout) :: this

    INTEGER :: no_children, nproc, i,j,iproc, maxsize
    INTEGER, ALLOCATABLE :: sizes(:,:)

    CHARACTER(len=*), PARAMETER :: routine = modname//':Cache_GPU_pointers'

    no_children = this%Get_no_of_children()
    IF (no_children == 0) CALL finish(routine, 'This should not happen')

    nproc = SIZE( this%mem )

    ALLOCATE(sizes(no_children,nproc))
    sizes(:,:) = 0

    !Note: process_action(iproc) = 2 means = AGGREGATE_ , but this can't be USEd
    !      from mo_jsb_process_class because of circular dependencies

    DO iproc=1,nproc
      IF (this%process_action(iproc) == 2) THEN
        IF ( ASSOCIATED(this%mem(iproc)%p) ) THEN
          DO i=1,no_children
            ! Collect the actual number of variables of all procs on all children
            IF ( ASSOCIATED(this%mem(iproc)%p%children(i)%p) ) &
              & sizes(i,iproc) = this%mem(iproc)%p%children(i)%p%no_of_vars
          END DO
        END IF
      END IF
    END DO
    maxsize = MAXVAL( sizes )

    ALLOCATE( this%ptrs2d_cache(no_children,maxsize,nproc), &
              this%ptrs3d_cache(no_children,maxsize,nproc)  )
    !$ACC ENTER DATA CREATE(this%ptrs2d_cache, this%ptrs3d_cache)

    DO iproc=1,nproc
      IF (this%process_action(iproc) == 2) THEN
        IF ( ASSOCIATED(this%mem(iproc)%p) ) THEN

          DO i=1,no_children
            IF ( ASSOCIATED(this%mem(iproc)%p%children(i)%p) ) THEN
              DO j=1,this%mem(iproc)%p%children(i)%p%no_of_vars
                IF ( ASSOCIATED(this%mem(iproc)%p%children(i)%p%vars(j)%p) ) THEN

                  IF (ASSOCIATED(this%mem(iproc)%p%children(i)%p%vars(j)%p%ptr2d)) THEN
                    IF ( .NOT. ASSOCIATED(this%ptrs2d_cache(i,j,iproc)%p, &
                                          this%mem(iproc)%p%children(i)%p%vars(j)%p%ptr2d) ) THEN
                      ! 2d cache pointer to proc iproc, child i, 2d variable j
                      this%ptrs2d_cache(i,j,iproc)%p => this%mem(iproc)%p%children(i)%p%vars(j)%p%ptr2d
                      __acc_attach(this%ptrs2d_cache(i,j,iproc)%p)
                    END IF
                  END IF

                  IF (ASSOCIATED(this%mem(iproc)%p%children(i)%p%vars(j)%p%ptr3d)) THEN
                    IF ( .NOT. ASSOCIATED(this%ptrs3d_cache(i,j,iproc)%p, &
                                          this%mem(iproc)%p%children(i)%p%vars(j)%p%ptr3d) ) THEN
                      this%ptrs3d_cache(i,j,iproc)%p => this%mem(iproc)%p%children(i)%p%vars(j)%p%ptr3d
                      __acc_attach(this%ptrs3d_cache(i,j,iproc)%p)
                    END IF
                  END IF

                END IF
              END DO
            END IF
          END DO
        END IF
      END IF
    END DO

  END SUBROUTINE Cache_GPU_pointers

#endif
END MODULE mo_jsb_tile_class
