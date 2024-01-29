! This module contains utilty types
!
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_grid_subset

  USE mo_kind,           ONLY: wp
  USE mo_exception,      ONLY: warning, finish, message
  USE mo_model_domain,   ONLY: t_patch, t_subset_range, t_subset_range_index, &
                               t_subset_indexed
  USE mo_mpi,            ONLY: get_my_mpi_work_id
  USE mo_impl_constants, ONLY: on_cells, on_edges, on_vertices
  USE mo_parallel_config, ONLY: nproma
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf

  IMPLICIT NONE

  PUBLIC :: t_subset_range, t_subset_range_index, t_subset_indexed

  PUBLIC :: fill_subset, get_index_range, fill_subset_indices
  PUBLIC :: read_subset, write_subset
  PUBLIC :: block_no, index_no, index_1d

CONTAINS

  !----------------------------------------------------
  !>
  ! Fills the subset_range with the indexes which statisfy
  !     start_mask <= mask <= end_mask
  !
  ! Assumes that invalid places in mask have value outside the
  ! start_mask - end_mask
  !PUBLIC ::
  ! Assumes that mask is of the shape (1:,1:)
  SUBROUTINE fill_subset(subset, patch, mask, start_mask, end_mask, subset_name, located)
    TYPE(t_subset_range), INTENT(inout) :: subset
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset
    INTEGER, OPTIONAL, INTENT(in) :: located
    CHARACTER(len=32), OPTIONAL :: subset_name
                                                ! unless the patch is declared INTENT(in)!
    INTEGER, INTENT(in) :: mask(:,:), start_mask, end_mask

    INTEGER :: masks_size(2)
    INTEGER :: block, index_in_block, start_index, end_index, start_block_loop

    CHARACTER(*), PARAMETER :: method_name = "mo_grid_subset:fill_subset"

    masks_size = SHAPE(mask)

    subset%start_block = -1
    subset%start_index = -1
    subset%end_block   = -2
    subset%end_index   = -2
    subset%block_size  = masks_size(1)
    subset%no_of_holes = 0
    subset%size        = 0
    subset%entity_location  = 0
    subset%patch       => patch
    subset%max_vertical_levels = 0
    NULLIFY(subset%vertical_levels)

    IF (PRESENT(located)) subset%entity_location = located

    DO block=1,masks_size(2)
      DO index_in_block=1,subset%block_size

        IF (mask(index_in_block, block) >= start_mask .AND. &
            mask(index_in_block, block) <= end_mask) THEN
          ! we found an element in range
          IF (subset%start_block < 0) THEN
            ! this is the first element
            subset%start_block = block
            subset%start_index = index_in_block
          ENDIF
          ! this is the up-to-now last element
          subset%end_block = block
          subset%end_index = index_in_block
        ENDIF

      ENDDO
    ENDDO

!     IF (subset%start_block == -1) THEN
!       CALL warning(method_name, "Empty range subset")
!     ENDIF

    IF (subset%start_block > -1) THEN
      ! count the holes
      DO block = subset%start_block, subset%end_block

        start_index = 1
        end_index = subset%block_size
        IF (block == subset%start_block) start_index = subset%start_index
        IF (block == subset%end_block)   end_index = subset%end_index

        DO index_in_block = start_index, end_index

          IF (mask(index_in_block, block) < start_mask .OR. &
              mask(index_in_block, block) > end_mask) THEN
            ! element is a hole in the range
            subset%no_of_holes = subset%no_of_holes + 1
          ENDIF

        ENDDO
      ENDDO
    ENDIF

    ! compute size
    subset%size = 0
    IF (subset%end_block > 0) THEN
      IF ((subset%end_block - subset%start_block) > 1) &
        subset%size = (subset%end_block - subset%start_block -1) * subset%block_size
      subset%size = subset%size + subset%end_index + (subset%block_size - subset%start_index + 1)
    ENDIF
    IF (PRESENT(subset_name)) &
      & subset%name = TRIM(subset_name)

    ! Guarantee somewhat safe extents if subset is empty. This should work in loops and array expressions.
    IF (subset%start_block < 0) THEN
      subset%start_block = 1
      subset%end_block = 0
      subset%start_index = 1
      subset%end_index = 0
    END IF
!     IF (subset%no_of_holes > 0) THEN
!       CALL warning(method_name, "We have holes in the range subset")
!     ENDIF

  END SUBROUTINE fill_subset

  !----------------------------------------------------
  ! fills the subset with indices according to the SR call
  ! (used to fill the subset owned_no_boundary which addresses
  !  prognostic cells only, excluding the boundary interpolation
  !  zone which is necessary for LAM or subdomain grids)
  !
  SUBROUTINE fill_subset_indices(subset, patch, start_block, end_block, start_index, end_index)
    TYPE(t_subset_range), INTENT(inout) :: subset
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset
    INTEGER, INTENT(in) :: start_block, end_block, start_index, end_index

    subset%start_block = start_block
    subset%end_block   = end_block
    subset%start_index = start_index
    subset%end_index   = end_index
    subset%block_size  = nproma
    subset%no_of_holes = -1
    subset%size        = -1
    subset%entity_location  = 0
    subset%patch       => patch
    subset%max_vertical_levels = 0

    NULLIFY(subset%vertical_levels)

  END SUBROUTINE fill_subset_indices
  !----------------------------------------------------

  !----------------------------------------------------

  SUBROUTINE get_index_range(subset_range, current_block, start_index, end_index)
    TYPE(t_subset_range), INTENT(in) :: subset_range
    INTEGER, INTENT(in) :: current_block
    INTEGER, INTENT(out) :: start_index, end_index
    !$ACC ROUTINE SEQ

    IF (current_block < subset_range%start_block .OR. &
        current_block > subset_range%end_block) THEN
      start_index = 1
      end_index = 0
    ELSE
      start_index = 1
      end_index = subset_range%block_size

      IF (current_block == subset_range%start_block) &
        start_index = subset_range%start_index
      IF (current_block == subset_range%end_block) &
        end_index = subset_range%end_index
    END IF

  END SUBROUTINE get_index_range
  !----------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE compare_subsets_range(subset_1, subset_2)
    TYPE(t_subset_range), INTENT(in) :: subset_1, subset_2

    CHARACTER(*), PARAMETER :: method_name = "mo_grid_subset:compare_subsets_range"

    IF (subset_1%entity_location /= subset_2%entity_location) THEN
      CALL message(method_name, TRIM(subset_1%name)//" "//TRIM(subset_1%name)//"subsets not of the same entity")
      RETURN
    ENDIF

  END SUBROUTINE compare_subsets_range
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! The following functions are for conversion of 1D to 2D indices and vice versa
  !
  ! Treatment of 0 (important for empty patches) and negative numbers:
  !
  ! Converting 1D => 2D:
  !
  ! 0 always is mapped to blk_no = 0, idx_no = 0
  ! negative numbers: Convert usings ABS(j) and negate idx_no
  !
  ! Thus: blk_no >= 0 always!
  !       idx_no > 0  for j > 0
  !       idx_no = 0  for j = 0
  !       idx_no < 0  for j < 0
  !
  ! This mimics mostly the behaviour of reshape_idx in mo_model_domimp_patches
  ! with a difference for nproma=1 and j=0 (where reshape_idx returns blk_no=0, idx_no=1)
  !
  ! The consisten treatment of 0 in the above way is very important for empty patches
  ! where start_index=1, end_index=0
  !
  !-------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION block_no(j)
    INTEGER, INTENT(in) :: j

    IF ( j == 0 ) THEN
      block_no = 0
    ELSE
      block_no = (ABS(j)-1)/nproma + 1
    ENDIF

  END FUNCTION block_no
  !-------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION index_no(j)
    INTEGER, INTENT(in) :: j

    IF( j == 0 ) THEN
      index_no = 0
    ELSE
      index_no = SIGN(MOD(ABS(j)-1,nproma)+1, j)
    ENDIF

  END FUNCTION index_no
  !-------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION index_1d(idx,block)
    INTEGER, INTENT(in) :: idx, block

    IF( block <= 0 ) THEN
      index_1d = 0 ! This covers the special case nproma==1,jb=0,jl=1
                 ! All other cases are invalid and get also a 0 returned
    ELSE
      index_1d = SIGN((block-1)*nproma + ABS(idx), idx)
    ENDIF

  END FUNCTION index_1d
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE read_subset(ncid, subset_range, patch)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_subset_range), INTENT(inout) :: subset_range
    TYPE(t_patch), TARGET, INTENT(in) :: patch  ! nag does not return the values in subset_range

    INTEGER :: netcd_status
    CHARACTER(*), PARAMETER :: method_name = "read_subset_range"

    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.start_block', subset_range%start_block)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read start_block")
    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.start_index', subset_range%start_index)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read start_index")

    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.end_block', subset_range%end_block)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read end_block")
    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.end_index', subset_range%end_index)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read end_index")

    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.block_size', subset_range%block_size)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read block_size")

    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.entity_location', subset_range%entity_location)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read entity_type")

    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.no_of_holes', subset_range%no_of_holes)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read no_of_holes")

    netcd_status = nfx_get_att(ncid, nf_global,TRIM(subset_range%name)//'.is_in_domain', subset_range%is_in_domain)
    IF (netcd_status /= nf_noerr) &
      & CALL finish(method_name, "Could not read is_in_domain")

    subset_range%patch => patch

  END SUBROUTINE read_subset
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE write_subset(ncid, subset_range)
    INTEGER, INTENT(in) :: ncid
    TYPE(t_subset_range), INTENT(inout) :: subset_range

    CHARACTER(*), PARAMETER :: method_name = "write_subset_range"

    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.start_block', nf_int,     &
      &  subset_range%start_block), method_name)
    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.start_index', nf_int,     &
      &  subset_range%start_index), method_name)

    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.end_block',   nf_int,     &
      & subset_range%end_block), method_name)
    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.end_index',   nf_int,     &
      & subset_range%end_index), method_name)

    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.block_size',  nf_int,     &
      & subset_range%block_size), method_name)

    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.entity_location', nf_int, &
      & subset_range%entity_location), method_name)

    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.no_of_holes', nf_int,     &
      & subset_range%no_of_holes), method_name)

    CALL nf(nfx_put_att(ncid, nf_global,TRIM(subset_range%name)//'.is_in_domain',nf_int,     &
      & subset_range%is_in_domain), method_name)

  END SUBROUTINE write_subset

END MODULE mo_grid_subset
