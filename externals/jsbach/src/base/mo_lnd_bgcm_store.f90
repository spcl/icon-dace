!> Summary: defines methods for collecting, storing and providing bgc materials
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
!>#### defines methods (used on tiles) for collecting, storing and providing bgc materials
!>
MODULE mo_lnd_bgcm_store
#ifndef __NO_QUINCY__

  USE mo_exception,             ONLY: finish, message
  USE mo_kind,                  ONLY: wp
  USE mo_jsb_control,           ONLY: debug_on
  USE mo_util,                  ONLY: int2string
  USE mo_lnd_bgcm_store_class,  ONLY: t_lnd_bgcm_store_abstract

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_lnd_bgcm_store, get_bgcm_idx


  ! ======================================================================================================= !
  !>
  !> bgc material store - stores pointers to the element variables of all bgcms across all processes
  !>
  TYPE, EXTENDS(t_lnd_bgcm_store_abstract) :: t_lnd_bgcm_store
  CONTAINS
    !> These deferred routines must match the interfaces provided in mo_lnd_bgcm_store_class
    PROCEDURE :: Init           => Init_bgcm_store
    PROCEDURE :: Append_bgcms   => Append_process_bgcms_to_store
    !> Additional procedures
    PROCEDURE :: Store_bgcms_in_matrices
    PROCEDURE :: Write_stored_matrices_to_bgcms
    PROCEDURE :: Get_matrix_from_store_2l_2d
    PROCEDURE :: Get_matrix_from_store_2l_3d
    PROCEDURE :: Get_matrix_from_store_1l_2d
    PROCEDURE :: Get_matrix_from_store_1l_3d

    FINAL   :: Finalize_lnd_bgcm_store
  END TYPE

  CHARACTER(len=*), PARAMETER :: modname = 'mo_lnd_bgcm_store'

CONTAINS

  ! ======================================================================================================= !
  !>
  !> Init this instance of a bgc material store
  !>
  SUBROUTINE Init_bgcm_store(this)
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store), INTENT(INOUT) :: this !< to be initialised bgcm store
    ! -------------------------------------------------------------------------------------------------- !

    ! After counting the materials the store can be prepared to
    ! - collect bookkeeping info about them
    ALLOCATE(this%bgc_material_name(this%nr_of_bgc_materials))
    ALLOCATE(this%bgc_material_ID(this%nr_of_bgc_materials))
    ALLOCATE(this%nr_of_levels(this%nr_of_bgc_materials))
    ALLOCATE(this%nr_of_dimensions(this%nr_of_bgc_materials))
    ALLOCATE(this%idx_in_store(this%nr_of_bgc_materials))

    ! - collect pointers to the bgc materials
    IF(this%nr_of_1l_2d_bgc_materials > 0) &
      & ALLOCATE(this%store_1l_2d_bgcms(this%nr_of_1l_2d_bgc_materials))
    IF(this%nr_of_1l_3d_bgc_materials > 0) &
      & ALLOCATE(this%store_1l_3d_bgcms(this%nr_of_1l_3d_bgc_materials))
    IF(this%nr_of_2l_2d_bgc_materials > 0) &
      & ALLOCATE(this%store_2l_2d_bgcms(this%nr_of_2l_2d_bgc_materials))
    IF(this%nr_of_2l_3d_bgc_materials > 0) &
      & ALLOCATE(this%store_2l_3d_bgcms(this%nr_of_2l_3d_bgc_materials))

  END SUBROUTINE Init_bgcm_store

  ! ======================================================================================================= !
  !>
  !> Append the bgcms of a process to this instance of a bgc material store
  !>
  SUBROUTINE Append_process_bgcms_to_store(this, main_bgcm, nblks, nproma, nsoil, tile_name, idx_bgcm, &
      &                                    idx_1l_2d_bgcm, idx_2l_2d_bgcm, idx_1l_3d_bgcm, idx_2l_3d_bgcm)

    USE mo_jsb_pool_class,      ONLY: t_jsb_pool
    USE mo_lnd_bgcm_class,      ONLY: get_dimension_of_element_variables
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store), INTENT(INOUT) :: this           !< this bgcm store instance
    CLASS(t_jsb_pool),       INTENT(IN)    :: main_bgcm      !< 'main' bgcm of a process
    CHARACTER(len=*),        INTENT(IN)    :: tile_name      !< name of the tile with this bgcm store
    INTEGER,                 INTENT(IN)    :: nblks          !< number of blocks
    INTEGER,                 INTENT(IN)    :: nproma         !< number of grid-cells
    INTEGER,                 INTENT(IN)    :: nsoil          !< number of soil layers
    INTEGER,                 INTENT(INOUT) :: idx_bgcm       !< idx of last added bgcm in the bgcm store
    INTEGER,                 INTENT(INOUT) :: idx_1l_2d_bgcm !< idx of last 1l 2d bgcm
    INTEGER,                 INTENT(INOUT) :: idx_2l_2d_bgcm !< idx of last 2l 2d bgcm
    INTEGER,                 INTENT(INOUT) :: idx_1l_3d_bgcm !< idx of last 1l 3d bgcm
    INTEGER,                 INTENT(INOUT) :: idx_2l_3d_bgcm !< idx of last 2l 3d bgcm
    ! -------------------------------------------------------------------------------------------------- !
    CLASS(t_jsb_pool),   POINTER  :: bgcm
    INTEGER                       :: i_bgcm, i_compartment, i_element
    INTEGER                       :: nr_of_bgcms_this_proc, nr_of_elements, nr_of_compartments

    CHARACTER(len=*), PARAMETER   :: routine = modname//':Append_process_bgcms_to_store'
    ! -------------------------------------------------------------------------------------------------- !

    IF (debug_on()) &
      & CALL message(TRIM(routine), 'Append bgcms: '//TRIM(main_bgcm%name)//' - tile '//TRIM(tile_name))

    nr_of_bgcms_this_proc = SIZE(main_bgcm%pool_list)
    DO i_bgcm = 1, nr_of_bgcms_this_proc
      idx_bgcm = idx_bgcm + 1
      bgcm => main_bgcm%pool_list(i_bgcm)%p

      this%bgc_material_name(idx_bgcm) = bgcm%name
      IF(bgcm%id == -1) THEN
        CALL finish(TRIM(routine), 'Found bgcm with no specified id '//TRIM(bgcm%name)// &
          & ' in tile '//TRIM(tile_name)//', please check!')
      ENDIF
      this%bgc_material_ID(idx_bgcm) = bgcm%id

      ! Currently only two cases are expected
      ! - either directly elements and no compartments (1l)
      ! - or compartments with elements but no direct elements (2l)
      IF(bgcm%contains_elements) THEN
        ! i.e. directly elements and no compartments (asserted in the counting routine)
        this%nr_of_levels(idx_bgcm) = 1
        this%nr_of_dimensions(idx_bgcm) = get_dimension_of_element_variables(bgcm, tile_name, routine)
        nr_of_elements = SIZE(bgcm%element_list)

        ! Currently two cases are expected: 2D or 3D variables
        SELECT CASE(this%nr_of_dimensions(idx_bgcm))
        CASE(2)
          idx_1l_2d_bgcm = idx_1l_2d_bgcm + 1
          ALLOCATE(this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%collection_1l_p(nr_of_elements))
          ALLOCATE(this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%mt_1l_2d_bgcm(nr_of_elements, nproma, nblks))
          this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%mt_1l_2d_bgcm = 0.0_wp
          this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%nr_of_elems = nr_of_elements
          this%idx_in_store(idx_bgcm) = idx_1l_2d_bgcm
          DO i_element = 1, nr_of_elements
            this%store_1l_2d_bgcms(idx_1l_2d_bgcm)%collection_1l_p(i_element)%p &
              & => bgcm%element_list(i_element)%p
          ENDDO
        CASE(3)
          idx_1l_3d_bgcm = idx_1l_3d_bgcm + 1
          ALLOCATE(this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%collection_1l_p(nr_of_elements))
          ALLOCATE(this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%mt_1l_3d_bgcm(nr_of_elements, nproma, nsoil, nblks))
          this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%mt_1l_3d_bgcm = 0.0_wp
          this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%nr_of_elems = nr_of_elements
          this%idx_in_store(idx_bgcm) = idx_1l_3d_bgcm
          DO i_element = 1, nr_of_elements
            this%store_1l_3d_bgcms(idx_1l_3d_bgcm)%collection_1l_p(i_element)%p &
              & => bgcm%element_list(i_element)%p
          ENDDO
        ENDSELECT

      ELSE
        ! i.e. no elements but one layer of compartments which are then expected to have elements
        ! but no further compartments (asserted in the counting routine)
        this%nr_of_levels(idx_bgcm) = 2
        ! current assumption (asserted in the counting routine): all elements of all compartments have the same dimension
        this%nr_of_dimensions(idx_bgcm) &
          & = get_dimension_of_element_variables(bgcm%pool_list(1)%p, tile_name, routine)

        nr_of_compartments = SIZE(bgcm%pool_list)
        nr_of_elements = SIZE(bgcm%pool_list(1)%p%element_list)

        SELECT CASE(this%nr_of_dimensions(idx_bgcm))
        CASE(2)
          idx_2l_2d_bgcm = idx_2l_2d_bgcm + 1
          ALLOCATE(this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%collection_2l_p(nr_of_compartments, nr_of_elements))
          ALLOCATE(this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%mt_2l_2d_bgcm(nr_of_compartments, nr_of_elements, &
            &                                                           nproma, nblks))
          this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%mt_2l_2d_bgcm = 0.0_wp
          this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%nr_of_elems = nr_of_elements
          this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%nr_of_parts = nr_of_compartments
          this%idx_in_store(idx_bgcm) = idx_2l_2d_bgcm
        CASE(3)
          idx_2l_3d_bgcm = idx_2l_3d_bgcm + 1
          ALLOCATE(this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%collection_2l_p(nr_of_compartments, nr_of_elements))
          ALLOCATE(this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%mt_2l_3d_bgcm(nr_of_compartments, nr_of_elements, &
            &                                                           nproma, nsoil, nblks))
          this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%mt_2l_3d_bgcm = 0.0_wp
          this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%nr_of_elems = nr_of_elements
          this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%nr_of_parts = nr_of_compartments
          this%idx_in_store(idx_bgcm) = idx_2l_3d_bgcm
        ENDSELECT

        DO i_compartment = 1, nr_of_compartments
          SELECT CASE(this%nr_of_dimensions(idx_bgcm))
          CASE(2)
            DO i_element = 1, nr_of_elements
              this%store_2l_2d_bgcms(idx_2l_2d_bgcm)%collection_2l_p(i_compartment, i_element)%p &
                & => bgcm%pool_list(i_compartment)%p%element_list(i_element)%p
            ENDDO
          CASE(3)
            DO i_element = 1, nr_of_elements
              this%store_2l_3d_bgcms(idx_2l_3d_bgcm)%collection_2l_p(i_compartment, i_element)%p &
                & => bgcm%pool_list(i_compartment)%p%element_list(i_element)%p
            ENDDO
          ENDSELECT
        ENDDO ! compartments
      ENDIF ! bgcm type
    ENDDO ! bgcms

  END SUBROUTINE Append_process_bgcms_to_store


  ! ======================================================================================================= !
  !>
  !> Store current values of all bgcms in the storage matrices
  !>
  SUBROUTINE Store_bgcms_in_matrices(this, tile_name, ics, ice, iblk)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store), INTENT(INOUT) :: this           !< this bgcm store instance
    CHARACTER(len=*),        INTENT(IN)    :: tile_name      !< name of the tile with this bgcm store
    INTEGER,                 INTENT(IN)    :: ics            !< start index
    INTEGER,                 INTENT(IN)    :: ice            !< end index
    INTEGER,                 INTENT(IN)    :: iblk           !< current block
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                :: i_bgcm         !< index of bgcm
    INTEGER                                :: i_part         !< index of compartment
    INTEGER                                :: i_elem         !< index of element

    CHARACTER(len=*), PARAMETER :: routine = modname//':Store_bgcms_in_matrices'
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on()) &
      & CALL message(TRIM(routine), 'Store bgcms in matrices for - tile '//TRIM(tile_name))

    DO i_bgcm = 1, this%nr_of_bgc_materials
      SELECT CASE(this%nr_of_levels(i_bgcm))
      CASE(1)
        ! i.e. only elements
        SELECT CASE(this%nr_of_dimensions(i_bgcm))
        CASE(2)
          ! 1level 2D -> fill matrix (elems, nc, nblks)
          DO i_elem = 1,this%store_1l_2d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
            this%store_1l_2d_bgcms(this%idx_in_store(i_bgcm))%mt_1l_2d_bgcm(i_elem,ics:ice,iblk) &
              & = this%store_1l_2d_bgcms(this%idx_in_store(i_bgcm))%collection_1l_p(i_elem)%p%ptr2d(ics:ice,iblk)
          ENDDO
        CASE(3)
          ! 1level 3D -> fill matrix (elems, nc, var_dim3, nblks)
          DO i_elem = 1,this%store_1l_3d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
            this%store_1l_3d_bgcms(this%idx_in_store(i_bgcm))%mt_1l_3d_bgcm(i_elem,ics:ice,:,iblk) &
              & = this%store_1l_3d_bgcms(this%idx_in_store(i_bgcm))%collection_1l_p(i_elem)%p%ptr3d(ics:ice,:,iblk)
          ENDDO
        ENDSELECT
      CASE(2)
        ! i.e. with compartments
        SELECT CASE(this%nr_of_dimensions(i_bgcm))
        CASE(2)
          ! 2levels 2D -> fill matrix (parts, elems, nc, nblks)
          DO i_part = 1,this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%nr_of_parts
            DO i_elem = 1,this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
              this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%mt_2l_2d_bgcm(i_part,i_elem,ics:ice,iblk) &
                & = this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%collection_2l_p(i_part, i_elem)%p%ptr2d(ics:ice,iblk)
            ENDDO
          ENDDO
        CASE(3)
          ! 2levels 3D -> fill matrix (parts, elems, nc, var_dim3, nblks)
          DO i_part = 1,this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%nr_of_parts
            DO i_elem = 1,this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
              this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%mt_2l_3d_bgcm(i_part,i_elem,ics:ice,:,iblk) &
                & = this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%collection_2l_p(i_part, i_elem)%p%ptr3d(ics:ice,:,iblk)
            ENDDO
          ENDDO
        ENDSELECT
      ENDSELECT
    ENDDO

  END SUBROUTINE Store_bgcms_in_matrices


  ! ======================================================================================================= !
  !>
  !> Write all stored matrices back to bgcms
  !>
  SUBROUTINE Write_stored_matrices_to_bgcms(this, tile_name, ics, ice, iblk)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store), INTENT(INOUT) :: this           !< this bgcm store instance
    CHARACTER(len=*),        INTENT(IN)    :: tile_name      !< name of the tile with this bgcm store
    INTEGER,                  INTENT(IN)   :: ics            !< start index
    INTEGER,                  INTENT(IN)   :: ice            !< end index
    INTEGER,                  INTENT(IN)   :: iblk           !< current block
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                :: i_bgcm         !< index of bgcm
    INTEGER                                :: i_part         !< index of compartment
    INTEGER                                :: i_elem         !< index of element

    CHARACTER(len=*), PARAMETER :: routine = modname//':Write_stored_matrices_to_bgcms'
    ! ----------------------------------------------------------------------------------------------------- !

    IF (debug_on()) &
      & CALL message(TRIM(routine), 'Store bgcms in matrices for - tile '//TRIM(tile_name))

    DO i_bgcm = 1, this%nr_of_bgc_materials
      SELECT CASE(this%nr_of_levels(i_bgcm))
      CASE(1)
        ! i.e. only elements
        SELECT CASE(this%nr_of_dimensions(i_bgcm))
        CASE(2)
          ! 1level 2D -> fill matrix (elems, nc, nblks)
          DO i_elem = 1,this%store_1l_2d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
            this%store_1l_2d_bgcms(this%idx_in_store(i_bgcm))%collection_1l_p(i_elem)%p%ptr2d(ics:ice,iblk) &
              & = this%store_1l_2d_bgcms(this%idx_in_store(i_bgcm))%mt_1l_2d_bgcm(i_elem,ics:ice,iblk)
          ENDDO
        CASE(3)
          ! 1level 3D -> fill matrix (elems, nc, var_dim3, nblks)
          DO i_elem = 1,this%store_1l_3d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
            this%store_1l_3d_bgcms(this%idx_in_store(i_bgcm))%collection_1l_p(i_elem)%p%ptr3d(ics:ice,:,iblk) &
              & = this%store_1l_3d_bgcms(this%idx_in_store(i_bgcm))%mt_1l_3d_bgcm(i_elem,ics:ice,:,iblk)
          ENDDO
        ENDSELECT
      CASE(2)
        ! i.e. with compartments
        SELECT CASE(this%nr_of_dimensions(i_bgcm))
        CASE(2)
          ! 2levels 2D -> fill matrix (parts, elems, nc, nblks)
          DO i_part = 1,this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%nr_of_parts
            DO i_elem = 1,this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
              this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%collection_2l_p(i_part, i_elem)%p%ptr2d(ics:ice,iblk) &
                & = this%store_2l_2d_bgcms(this%idx_in_store(i_bgcm))%mt_2l_2d_bgcm(i_part,i_elem,ics:ice,iblk)
            ENDDO
          ENDDO
        CASE(3)
          ! 2levels 3D -> fill matrix (parts, elems, nc, var_dim3, nblks)
          DO i_part = 1,this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%nr_of_parts
            DO i_elem = 1,this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%nr_of_elems
              this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%collection_2l_p(i_part, i_elem)%p%ptr3d(ics:ice,:,iblk) &
                & = this%store_2l_3d_bgcms(this%idx_in_store(i_bgcm))%mt_2l_3d_bgcm(i_part,i_elem,ics:ice,:,iblk)
            ENDDO
          ENDDO
        ENDSELECT
      ENDSELECT
    ENDDO

  END SUBROUTINE Write_stored_matrices_to_bgcms

  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) from the store for this 2l 2d bgcm
  !>
  FUNCTION Get_matrix_from_store_2l_2d(this, bgcm_id, ics, ice, iblk, tile_name) RESULT(bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this               !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id            !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                !< start index
    INTEGER,                    INTENT(IN)    :: ice                !< end index
    INTEGER,                    INTENT(IN)    :: iblk               !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name          !< name of the current tile
    REAL(wp),                   POINTER       :: bgcm_matrix(:,:,:) !< dims: (part,elem,nc)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_store_2l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Get index of the queried bgcm in this store
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! 2levels 2D -> fill matrix (parts, elems, nc)
    bgcm_matrix => this%store_2l_2d_bgcms(this%idx_in_store(bgcm_idx))%mt_2l_2d_bgcm(:,:,ics:ice,iblk)

  END FUNCTION Get_matrix_from_store_2l_2d


  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) from the store for this 2l 3d bgcm
  !>
  FUNCTION Get_matrix_from_store_2l_3d(this, bgcm_id, ics, ice, iblk, tile_name) RESULT(bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this                 !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id              !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                  !< start index
    INTEGER,                    INTENT(IN)    :: ice                  !< end index
    INTEGER,                    INTENT(IN)    :: iblk                 !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name            !< name of the current tile
    REAL(wp),                   POINTER       :: bgcm_matrix(:,:,:,:) !< dims: (part,elem,nc,layer)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_store_2l_3d'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Get index of the queried bgcm in this store
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! 2levels 3D -> fill matrix (parts, elems, nc, var_dim3)
    bgcm_matrix => this%store_2l_3d_bgcms(this%idx_in_store(bgcm_idx))%mt_2l_3d_bgcm(:, :,ics:ice,:,iblk)

  END FUNCTION Get_matrix_from_store_2l_3d


  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) from the store for this 1l 2d bgcm
  !>
  FUNCTION Get_matrix_from_store_1l_2d(this, bgcm_id, ics, ice, iblk, tile_name) RESULT(bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this             !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id          !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics              !< start index
    INTEGER,                    INTENT(IN)    :: ice              !< end index
    INTEGER,                    INTENT(IN)    :: iblk             !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name        !< name of the current tile
    REAL(wp),                   POINTER       :: bgcm_matrix(:,:) !< dims: (elem,nc)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_store_1l_2d'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Get index of the queried bgcm in this store
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! 1level 2D -> fill matrix (elems, nc)
    bgcm_matrix => this%store_1l_2d_bgcms(this%idx_in_store(bgcm_idx))%mt_1l_2d_bgcm(:,ics:ice,iblk)

  END FUNCTION Get_matrix_from_store_1l_2d


  ! ======================================================================================================= !
  !>
  !> Get a matrix for the chunk (ics:ice, iblk) from the store for this 1l 3d bgcm
  !>
  FUNCTION Get_matrix_from_store_1l_3d(this, bgcm_id, ics, ice, iblk, tile_name) RESULT(bgcm_matrix)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(INOUT) :: this               !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id            !< id of bgcm queried from the store
    INTEGER,                    INTENT(IN)    :: ics                !< start index
    INTEGER,                    INTENT(IN)    :: ice                !< end index
    INTEGER,                    INTENT(IN)    :: iblk               !< current block
    CHARACTER (len=*),          INTENT(IN)    :: tile_name          !< name of the current tile
    REAL(wp),                   POINTER       :: bgcm_matrix(:,:,:) !< dims: (elem,nc,layer)
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER                                   :: bgcm_idx           !< index of bgcm queried from the store

    CHARACTER(len=*), PARAMETER :: routine = modname//':Get_matrix_from_store_1l_3d'
    ! ----------------------------------------------------------------------------------------------------- !
    ! Get index of the queried bgcm in this store
    bgcm_idx = get_bgcm_idx(this, bgcm_id, tile_name, routine)

    ! 1level 3D -> fill matrix (elems, nc, var_dim3)
    bgcm_matrix => this%store_1l_3d_bgcms(this%idx_in_store(bgcm_idx))%mt_1l_3d_bgcm(:,ics:ice,:,iblk)

  END FUNCTION Get_matrix_from_store_1l_3d


  ! ======================================================================================================= !
  !>
  !> Get the index of this bgcm (id) from the store
  !>
  FUNCTION get_bgcm_idx(this, bgcm_id, tile_name, caller) RESULT(bgcm_idx)
    ! ----------------------------------------------------------------------------------------------------- !
    CLASS(t_lnd_bgcm_store),    INTENT(IN)    :: this        !< this bgcm store instance
    INTEGER,                    INTENT(IN)    :: bgcm_id     !< id of bgcm queried from the store
    CHARACTER (len=*),          INTENT(IN)    :: tile_name   !< name of the current tile
    CHARACTER(len=*),           INTENT(in)    :: caller      !< calling routine
    INTEGER                                   :: bgcm_idx    !< index of the queried bgcm in the store
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: idx       !< current index
    ! ----------------------------------------------------------------------------------------------------- !
    ! Assert a valid store
    IF(.NOT. ALLOCATED(this%bgc_material_ID)) CALL finish(TRIM(caller), 'Attempt to get bgcm from an invalid store!')

    bgcm_idx = 0
    DO idx = 1,SIZE(this%bgc_material_ID)
      IF (this%bgc_material_ID(idx) == bgcm_id) THEN
        IF (bgcm_idx /= 0) THEN
          CALL finish(TRIM(caller), 'Found not unique id '//int2string(bgcm_id)//      &
            & ' used for '//TRIM(this%bgc_material_name(bgcm_idx))//' and '//          &
            & TRIM(this%bgc_material_name(idx))//                                      &
            & ' in bgcm store of tile '//TRIM(tile_name)//', please check!')
        ENDIF
        bgcm_idx = idx
      ENDIF
    ENDDO

    IF (bgcm_idx == 0) THEN
      CALL finish(TRIM(caller), 'Did not find bgcm with given id '//int2string(bgcm_id)// &
        & ' in bgcm store of tile '//TRIM(tile_name)//', please check!')
    ENDIF

  END FUNCTION get_bgcm_idx


  ! ======================================================================================================= !
  !>
  !> deallocate bgcm store components
  !>
  SUBROUTINE Finalize_lnd_bgcm_store(bgcm_store)
    ! ----------------------------------------------------------------------------------------------------- !
    TYPE(t_lnd_bgcm_store) :: bgcm_store
    ! ----------------------------------------------------------------------------------------------------- !
    INTEGER :: i
    ! ----------------------------------------------------------------------------------------------------- !

    IF(ALLOCATED(bgcm_store%bgc_material_ID)) DEALLOCATE(bgcm_store%bgc_material_ID)
    IF(ALLOCATED(bgcm_store%nr_of_levels)) DEALLOCATE(bgcm_store%nr_of_levels)
    IF(ALLOCATED(bgcm_store%nr_of_dimensions)) DEALLOCATE(bgcm_store%nr_of_dimensions)
    IF(ALLOCATED(bgcm_store%idx_in_store)) DEALLOCATE(bgcm_store%idx_in_store)

    IF(ALLOCATED(bgcm_store%store_1l_2d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_1l_2d_bgcms)
        IF(ALLOCATED(bgcm_store%store_1l_2d_bgcms(i)%collection_1l_p)) &
          & DEALLOCATE(bgcm_store%store_1l_2d_bgcms(i)%collection_1l_p)
        IF(ASSOCIATED(bgcm_store%store_1l_2d_bgcms(i)%mt_1l_2d_bgcm)) &
          & DEALLOCATE(bgcm_store%store_1l_2d_bgcms(i)%mt_1l_2d_bgcm)
        IF(ASSOCIATED(bgcm_store%store_1l_2d_bgcms(i)%mt_1l_3d_bgcm)) &
          & DEALLOCATE(bgcm_store%store_1l_2d_bgcms(i)%mt_1l_3d_bgcm)
      ENDDO
      DEALLOCATE(bgcm_store%store_1l_2d_bgcms)
    ENDIF

    IF(ALLOCATED(bgcm_store%store_1l_3d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_1l_3d_bgcms)
        IF(ALLOCATED(bgcm_store%store_1l_3d_bgcms(i)%collection_1l_p)) &
          & DEALLOCATE(bgcm_store%store_1l_3d_bgcms(i)%collection_1l_p)
        IF(ASSOCIATED(bgcm_store%store_1l_3d_bgcms(i)%mt_1l_2d_bgcm)) &
          & DEALLOCATE(bgcm_store%store_1l_3d_bgcms(i)%mt_1l_2d_bgcm)
        IF(ASSOCIATED(bgcm_store%store_1l_3d_bgcms(i)%mt_1l_3d_bgcm)) &
          & DEALLOCATE(bgcm_store%store_1l_3d_bgcms(i)%mt_1l_3d_bgcm)
      ENDDO
      DEALLOCATE(bgcm_store%store_1l_3d_bgcms)
    ENDIF

    IF(ALLOCATED(bgcm_store%store_2l_2d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_2l_2d_bgcms)
        IF(ALLOCATED(bgcm_store%store_2l_2d_bgcms(i)%collection_2l_p)) &
          & DEALLOCATE(bgcm_store%store_2l_2d_bgcms(i)%collection_2l_p)
         IF(ASSOCIATED(bgcm_store%store_2l_2d_bgcms(i)%mt_2l_2d_bgcm)) &
           & DEALLOCATE(bgcm_store%store_2l_2d_bgcms(i)%mt_2l_2d_bgcm)
         IF(ASSOCIATED(bgcm_store%store_2l_2d_bgcms(i)%mt_2l_3d_bgcm)) &
           & DEALLOCATE(bgcm_store%store_2l_2d_bgcms(i)%mt_2l_3d_bgcm)
      ENDDO
      DEALLOCATE(bgcm_store%store_2l_2d_bgcms)
    ENDIF

    IF(ALLOCATED(bgcm_store%store_2l_3d_bgcms)) THEN
      DO i = 1,SIZE(bgcm_store%store_2l_3d_bgcms)
        IF(ALLOCATED(bgcm_store%store_2l_3d_bgcms(i)%collection_2l_p)) &
          & DEALLOCATE(bgcm_store%store_2l_3d_bgcms(i)%collection_2l_p)
        IF(ASSOCIATED(bgcm_store%store_2l_3d_bgcms(i)%mt_2l_2d_bgcm)) &
          & DEALLOCATE(bgcm_store%store_2l_3d_bgcms(i)%mt_2l_2d_bgcm)
        IF(ASSOCIATED(bgcm_store%store_2l_3d_bgcms(i)%mt_2l_3d_bgcm)) &
          & DEALLOCATE(bgcm_store%store_2l_3d_bgcms(i)%mt_2l_3d_bgcm)
      ENDDO
      DEALLOCATE(bgcm_store%store_2l_3d_bgcms)
    ENDIF

  END SUBROUTINE Finalize_lnd_bgcm_store


#endif
END MODULE mo_lnd_bgcm_store
