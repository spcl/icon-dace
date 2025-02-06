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

! Utility module: reorder cells, edges, vertices of a patch according
! to a given permutation array.

MODULE mo_reorder_patches

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: finish
  USE mo_parallel_config, ONLY: nproma
  USE mo_communication,   ONLY: idx_no, blk_no, idx_1d
  USE mo_model_domain,    ONLY: t_patch, t_tangent_vectors
  USE mo_math_types,      ONLY: t_geographical_coordinates, t_cartesian_coordinates
  USE mo_decomposition_tools, ONLY: t_grid_domain_decomp_info, &
    &                               t_glb2loc_index_lookup
  USE mo_util_sort,       ONLY: quicksort

  IMPLICIT NONE
  PRIVATE


  CHARACTER(LEN=*), PARAMETER :: modname = TRIM('mo_reorder_patches')

  PUBLIC :: reorder_cells
  PUBLIC :: reorder_edges
  PUBLIC :: reorder_verts

  INTERFACE reorder_array_pos
    MODULE PROCEDURE reorder_array_pos_i1D
    MODULE PROCEDURE reorder_array_pos_i2D
    MODULE PROCEDURE reorder_array_pos_l2D
    MODULE PROCEDURE reorder_array_pos_r2D
    MODULE PROCEDURE reorder_array_pos_g2D
    MODULE PROCEDURE reorder_array_pos_t2D
    MODULE PROCEDURE reorder_array_pos_c2D
    MODULE PROCEDURE reorder_array_pos_i3D
    MODULE PROCEDURE reorder_array_pos_r3D
    MODULE PROCEDURE reorder_array_pos_t3D
    MODULE PROCEDURE reorder_array_pos_i4D
  END INTERFACE

  INTERFACE reorder_array_content
    MODULE PROCEDURE reorder_array_content_i1D
    MODULE PROCEDURE reorder_array_content_i2D
    MODULE PROCEDURE reorder_array_content_i3D
    MODULE PROCEDURE reorder_array_content_i4D
  END INTERFACE

CONTAINS

  !> In-situ reordering of patch CELLS according to a given
  !  permutation array.
  !
  !  @note The permutation is applied to the local patch subset only,
  !        this subroutine is not MPI-collective!
  !
  !  @note we do not reorder the start_idx/blk and end_idx/blk data!
  !
  SUBROUTINE reorder_cells(pp, idx_old2new, opt_parent_pp, opt_child_pp)
    TYPE(t_patch),         INTENT(INOUT) :: pp               ! patch data structure
    INTEGER,               INTENT(IN)    :: idx_old2new(:)   ! permutation array
    TYPE(t_patch), INTENT(INOUT), OPTIONAL :: opt_parent_pp  ! (optional:) parent patch data structure
    TYPE(t_patch), INTENT(INOUT), OPTIONAL :: opt_child_pp   ! (optional:) child patch data structure
    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_cells"

    ! ----------------------------------------
    ! in this patch: translate array positions
    ! ----------------------------------------

    CALL reorder_array_pos(pp%cells%num_edges,       idx_old2new,      pp%nblks_c, pp%npromz_c)
    IF (ALLOCATED(pp%cells%parent_glb_idx)) &
    CALL reorder_array_pos(pp%cells%parent_glb_idx,  idx_old2new,      pp%nblks_c, pp%npromz_c)
    IF (ALLOCATED(pp%cells%parent_glb_blk)) &
    CALL reorder_array_pos(pp%cells%parent_glb_blk,  idx_old2new,      pp%nblks_c, pp%npromz_c)
    IF (ALLOCATED(pp%cells%child_idx)) &
    CALL reorder_array_pos(pp%cells%child_idx,       idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    IF (ALLOCATED(pp%cells%child_blk)) &
    CALL reorder_array_pos(pp%cells%child_blk,       idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%child_id,        idx_old2new,      pp%nblks_c, pp%npromz_c)
    CALL reorder_array_pos(pp%cells%phys_id,         idx_old2new,      pp%nblks_c, pp%npromz_c)
    CALL reorder_array_pos(pp%cells%edge_idx,        idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%edge_blk,        idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%vertex_idx,      idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%vertex_blk,      idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%edge_orientation,idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%refin_ctrl,      idx_old2new,      pp%nblks_c, pp%npromz_c)
    IF (ALLOCATED(pp%cells%ddqz_z_full)) THEN
      CALL reorder_array_pos(pp%cells%ddqz_z_full,     idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 3)
    END IF
    CALL reorder_array_pos(pp%cells%f_c,             idx_old2new,      pp%nblks_c, pp%npromz_c)
    CALL reorder_array_pos(pp%cells%area,            idx_old2new,      pp%nblks_c, pp%npromz_c)
    CALL reorder_array_pos(pp%cells%center,          idx_old2new,      pp%nblks_c, pp%npromz_c)
    CALL reorder_array_pos(pp%cells%cartesian_center,idx_old2new,      pp%nblks_c, pp%npromz_c)
    CALL reorder_decomp_info(pp%cells%decomp_info, idx_old2new, &
      &                      pp%n_patch_cells, pp%nblks_c, pp%npromz_c)

    ! ----------------------------------------------
    ! in this patch: translate position and contents
    ! ----------------------------------------------

    CALL reorder_array_content(pp%cells%neighbor_idx, pp%cells%neighbor_blk, idx_old2new, &
      &                        pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%neighbor_idx,    idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)
    CALL reorder_array_pos(pp%cells%neighbor_blk,    idx_old2new,      pp%nblks_c, pp%npromz_c, 1, 2)

    ! in this patch: translate contents of edges data structure
    CALL reorder_array_content(pp%edges%cell_idx, pp%edges%cell_blk, idx_old2new, &
      &                        pp%nblks_e, pp%npromz_e, 1, 2)
    CALL reorder_array_content(pp%edges%butterfly_idx, pp%edges%butterfly_blk, idx_old2new, &
      &                        pp%nblks_e, pp%npromz_e, 1, 2, opt_lcatch_zeros=.TRUE.)
    ! in this patch: translate contents of verts data structure
    CALL reorder_array_content(pp%verts%cell_idx, pp%verts%cell_blk, idx_old2new, &
      &                        pp%nblks_v, pp%npromz_v, 1, 2)

    ! -----------------------------------
    ! in parent patch: translate contents
    ! -----------------------------------

    IF (PRESENT(opt_parent_pp)) THEN
      CALL reorder_array_content(opt_parent_pp%cells%child_idx, opt_parent_pp%cells%child_blk, &
        &                        idx_old2new, opt_parent_pp%nblks_c, opt_parent_pp%npromz_c, 1, 2)
    END IF

    ! ----------------------------------
    ! in child patch: translate contents
    ! ----------------------------------

    IF (PRESENT(opt_child_pp)) THEN
      CALL reorder_array_content(opt_child_pp%cells%parent_loc_idx, &
        &                        opt_child_pp%cells%parent_loc_blk, &
        &                        idx_old2new, opt_child_pp%nblks_c, &
        &                        opt_child_pp%npromz_c)
    END IF

  END SUBROUTINE reorder_cells


  !> In-situ reordering of patch EDGES according to a given
  !  permutation array.
  !
  !  @note The permutation is applied to the local patch subset only,
  !        this subroutine is not MPI-collective!
  !
  SUBROUTINE reorder_edges(pp, idx_old2new, opt_parent_pp, opt_child_pp)
    TYPE(t_patch),         INTENT(INOUT) :: pp             ! patch data structure
    INTEGER,               INTENT(IN)    :: idx_old2new(:) ! permutation array
    TYPE(t_patch), INTENT(INOUT), OPTIONAL :: opt_parent_pp  ! (optional:) parent patch data structure
    TYPE(t_patch), INTENT(INOUT), OPTIONAL :: opt_child_pp   ! (optional:) child patch data structure

    ! ----------------------------------------
    ! in this patch: translate array positions
    ! ----------------------------------------
    IF (ALLOCATED(pp%edges%parent_glb_idx)) &
    CALL reorder_array_pos(pp%edges%parent_glb_idx,         idx_old2new,      pp%nblks_e, pp%npromz_e)
    IF (ALLOCATED(pp%edges%parent_glb_blk)) &
    CALL reorder_array_pos(pp%edges%parent_glb_blk,         idx_old2new,      pp%nblks_e, pp%npromz_e)
    IF (ALLOCATED(pp%edges%child_idx)) &
    CALL reorder_array_pos(pp%edges%child_idx,              idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    IF (ALLOCATED(pp%edges%child_blk)) &
    CALL reorder_array_pos(pp%edges%child_blk,              idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%child_id,               idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%phys_id,                idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%cell_idx,               idx_old2new,      pp%nblks_e, pp%npromz_e, 1, 2)
    CALL reorder_array_pos(pp%edges%cell_blk,               idx_old2new,      pp%nblks_e, pp%npromz_e, 1, 2)
    CALL reorder_array_pos(pp%edges%vertex_idx,             idx_old2new,      pp%nblks_e, pp%npromz_e, 1, 2)
    CALL reorder_array_pos(pp%edges%vertex_blk,             idx_old2new,      pp%nblks_e, pp%npromz_e, 1, 2)
    CALL reorder_array_pos(pp%edges%tangent_orientation,     idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%quad_idx,               idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%quad_blk,               idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%quad_orientation,       idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%butterfly_idx,          idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%butterfly_blk,          idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%center,                 idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%primal_normal,          idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%primal_cart_normal,     idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%dual_cart_normal,       idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%dual_normal,            idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%primal_normal_cell,     idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%dual_normal_cell,       idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%primal_normal_vert,     idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%dual_normal_vert,       idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%primal_edge_length,     idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%inv_primal_edge_length, idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%dual_edge_length,       idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%inv_dual_edge_length,   idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%edge_vert_length,       idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%edge_cell_length,       idx_old2new,      pp%nblks_e, pp%npromz_e, 1,2)
    CALL reorder_array_pos(pp%edges%inv_vert_vert_length,   idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%area_edge,              idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%quad_area,              idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%cartesian_center,       idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%cartesian_dual_middle,  idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%f_e,                    idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%fn_e,                   idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%ft_e,                   idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_array_pos(pp%edges%refin_ctrl,             idx_old2new,      pp%nblks_e, pp%npromz_e)
    CALL reorder_decomp_info(pp%edges%decomp_info, idx_old2new, &
      &                      pp%n_patch_edges, pp%nblks_e, pp%npromz_e)

    ! ----------------------------------------------
    ! in this patch: translate contents
    ! ----------------------------------------------

    CALL reorder_array_content(pp%edges%quad_idx, pp%edges%quad_blk, idx_old2new, &
      &                        pp%nblks_e, pp%npromz_e, 1, 2, opt_lcatch_zeros=.TRUE.)
    ! in this patch: translate contents of cells data structure
    CALL reorder_array_content(pp%cells%edge_idx, pp%cells%edge_blk, idx_old2new, &
      &                        pp%nblks_c, pp%npromz_c, 1, 2)
    ! in this patch: translate contents of verts data structure
    CALL reorder_array_content(pp%verts%edge_idx, pp%verts%edge_blk, idx_old2new, &
      &                        pp%nblks_v, pp%npromz_v, 1, 2)

    ! -----------------------------------
    ! in parent patch: translate contents
    ! -----------------------------------

    IF (PRESENT(opt_parent_pp)) THEN
      CALL reorder_array_content(opt_parent_pp%edges%child_idx, opt_parent_pp%edges%child_blk, &
        &                        idx_old2new, opt_parent_pp%nblks_e, opt_parent_pp%npromz_e, 1, 2)
    END IF

    ! ----------------------------------
    ! in child patch: translate contents
    ! ----------------------------------

    IF (PRESENT(opt_child_pp)) THEN
      CALL reorder_array_content(opt_child_pp%edges%parent_loc_idx, &
        &                        opt_child_pp%edges%parent_loc_blk, &
        &                        idx_old2new, opt_child_pp%nblks_e, &
        &                        opt_child_pp%npromz_e)
    END IF

  END SUBROUTINE reorder_edges


  !> In-situ reordering of patch VERTICES according to a given
  !  permutation array.
  !
  !  @note The permutation is applied to the local patch subset only,
  !        this subroutine is not MPI-collective!
  !
  SUBROUTINE reorder_verts(pp, idx_old2new)
    TYPE(t_patch),         INTENT(INOUT) :: pp             ! patch data structure
    INTEGER,               INTENT(IN)    :: idx_old2new(:) ! permutation array

    ! ----------------------------------------
    ! in this patch: translate array positions
    ! ----------------------------------------

    CALL reorder_array_pos(pp%verts%phys_id,                idx_old2new,      pp%nblks_v, pp%npromz_v)
    CALL reorder_array_pos(pp%verts%cell_idx,               idx_old2new,      pp%nblks_v, pp%npromz_v, 1,2)
    CALL reorder_array_pos(pp%verts%cell_blk,               idx_old2new,      pp%nblks_v, pp%npromz_v, 1,2)
    CALL reorder_array_pos(pp%verts%edge_idx,               idx_old2new,      pp%nblks_v, pp%npromz_v, 1,2)
    CALL reorder_array_pos(pp%verts%edge_blk,               idx_old2new,      pp%nblks_v, pp%npromz_v, 1,2)
    CALL reorder_array_pos(pp%verts%edge_orientation,       idx_old2new,      pp%nblks_v, pp%npromz_v, 1,2)
    CALL reorder_array_pos(pp%verts%num_edges,              idx_old2new,      pp%nblks_v, pp%npromz_v)
    CALL reorder_array_pos(pp%verts%vertex,                 idx_old2new,      pp%nblks_v, pp%npromz_v)
    CALL reorder_array_pos(pp%verts%dual_area,              idx_old2new,      pp%nblks_v, pp%npromz_v)
    CALL reorder_array_pos(pp%verts%f_v,                    idx_old2new,      pp%nblks_v, pp%npromz_v)
    CALL reorder_array_pos(pp%verts%cartesian,              idx_old2new,      pp%nblks_v, pp%npromz_v)
    CALL reorder_array_pos(pp%verts%refin_ctrl,             idx_old2new,      pp%nblks_v, pp%npromz_v)
    CALL reorder_decomp_info(pp%verts%decomp_info, idx_old2new, &
      &                      pp%n_patch_verts, pp%nblks_v, pp%npromz_v)

    ! ----------------------------------------------
    ! in this patch: translate position and contents
    ! ----------------------------------------------

    CALL reorder_array_content(pp%verts%neighbor_idx, pp%verts%neighbor_blk, idx_old2new, &
      &                        pp%nblks_v, pp%npromz_v, 1, 2, opt_lcatch_zeros=.TRUE.)
    ! in this patch: translate contents of cells data structure
    CALL reorder_array_content(pp%cells%vertex_idx, pp%cells%vertex_blk, idx_old2new,     &
      &                        pp%nblks_c, pp%npromz_c, 1, 2)
    ! in this patch: translate contents of edges data structure
    CALL reorder_array_content(pp%edges%vertex_idx, pp%edges%vertex_blk, idx_old2new,     &
      &                        pp%nblks_e, pp%npromz_e, 1, 2, opt_lcatch_zeros=.TRUE.)

  END SUBROUTINE reorder_verts

  SUBROUTINE reorder_decomp_info(decomp_info, idx_old2new, n, nblks, npromz)

    TYPE(t_grid_domain_decomp_info), INTENT(inout) :: decomp_info
    INTEGER, INTENT(IN), CONTIGUOUS :: idx_old2new(:) ! permutation array
    INTEGER, INTENT(IN) :: n, nblks, npromz

    CALL reorder_array_pos(decomp_info%owner_mask, idx_old2new, nblks, npromz)
    CALL reorder_array_pos(decomp_info%decomp_domain, idx_old2new, nblks, npromz)
    CALL reorder_array_pos(decomp_info%owner_local, idx_old2new, n)
    CALL reorder_array_pos(decomp_info%glb_index, idx_old2new, n)

    CALL reorder_glb2loc_index_lookup(decomp_info%glb2loc_index, &
      &                               decomp_info%glb_index, n)

  END SUBROUTINE reorder_decomp_info

  SUBROUTINE reorder_glb2loc_index_lookup(glb2loc_index, glb_index, n)

    TYPE(t_glb2loc_index_lookup), INTENT(inout) :: glb2loc_index
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: glb_index(n)

    INTEGER :: i

    DEALLOCATE(glb2loc_index%inner_glb_index, &
      &        glb2loc_index%inner_glb_index_to_loc, &
      &        glb2loc_index%outer_glb_index, &
      &        glb2loc_index%outer_glb_index_to_loc)
    ALLOCATE(glb2loc_index%inner_glb_index(0), &
      &      glb2loc_index%inner_glb_index_to_loc(0), &
      &      glb2loc_index%outer_glb_index(n), &
      &      glb2loc_index%outer_glb_index_to_loc(n))

    glb2loc_index%outer_glb_index(:) = glb_index(:)
    glb2loc_index%outer_glb_index_to_loc(:) = (/(i, i = 1, n)/)

    CALL quicksort(glb2loc_index%outer_glb_index(:), &
      &            glb2loc_index%outer_glb_index_to_loc(:))

  END SUBROUTINE reorder_glb2loc_index_lookup

  !> reorder array entries according to a given permutation array.
  !  2D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_i1D(arr, idx_old2new, nentries)
    INTEGER,     INTENT(INOUT), CONTIGUOUS :: arr(:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:) ! permutation array
    INTEGER,     INTENT(IN)    :: nentries

    ! local variables
    INTEGER :: tmp(SIZE(arr,1))
    INTEGER :: j

    tmp(:) = arr(:)
    DO j=1,nentries
      arr(idx_old2new(j)) = tmp(j)
    END DO
  END SUBROUTINE reorder_array_pos_i1D


  !> reorder array entries according to a given permutation array.
  !  2D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_i2D(arr, idx_old2new, nblks, npromz)
    INTEGER,     INTENT(INOUT), CONTIGUOUS :: arr(:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:) ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz

    ! local variables
    INTEGER :: tmp(SIZE(arr,1), SIZE(arr,2))
    INTEGER :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:) = arr(:,:)
    DO jb=1,nblks
      i_endidx = nproma
      if (jb == nblks) i_endidx = npromz
      DO jc=1,i_endidx
        iidx = idx_old2new(idx_1d(jc,jb))
        jc_new = idx_no(iidx)
        jb_new = blk_no(iidx)
        arr(jc_new,jb_new) = tmp(jc,jb)
      END DO
    END DO
  END SUBROUTINE reorder_array_pos_i2D


  !> reorder array entries according to a given permutation array.
  !  2D implementation, LOGICALS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_l2D(arr, idx_old2new, nblks, npromz)
    LOGICAL,     INTENT(INOUT), CONTIGUOUS :: arr(:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:) ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz

    ! local variables
    LOGICAL :: tmp(SIZE(arr,1), SIZE(arr,2))
    INTEGER :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:) = arr(:,:)
    DO jb=1,nblks
      i_endidx = nproma
      if (jb == nblks) i_endidx = npromz
      DO jc=1,i_endidx
        iidx = idx_old2new(idx_1d(jc,jb))
        jc_new = idx_no(iidx)
        jb_new = blk_no(iidx)
        arr(jc_new,jb_new) = tmp(jc,jb)
      END DO
    END DO
  END SUBROUTINE reorder_array_pos_l2D


  !> reorder array entries according to a given permutation array.
  !  2D implementation, REALS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_r2D(arr, idx_old2new, nblks, npromz)
    REAL(wp),    INTENT(INOUT), CONTIGUOUS :: arr(:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:) ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz

    ! local variables
    REAL(wp) :: tmp(SIZE(arr,1), SIZE(arr,2))
    INTEGER  :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:) = arr(:,:)
    DO jb=1,nblks
      i_endidx = nproma
      if (jb == nblks) i_endidx = npromz
      DO jc=1,i_endidx
        iidx = idx_old2new(idx_1d(jc,jb))
        jc_new = idx_no(iidx)
        jb_new = blk_no(iidx)
        arr(jc_new,jb_new) = tmp(jc,jb)
      END DO
    END DO
  END SUBROUTINE reorder_array_pos_r2D


  !> reorder array entries according to a given permutation array.
  !  2D implementation, geographical coordinates.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_g2D(arr, idx_old2new, nblks, npromz)
    TYPE(t_geographical_coordinates), INTENT(INOUT), CONTIGUOUS :: arr(:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:) ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz

    ! local variables
    TYPE(t_geographical_coordinates) :: tmp(SIZE(arr,1), SIZE(arr,2))
    INTEGER  :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:) = arr(:,:)
    DO jb=1,nblks
      i_endidx = nproma
      if (jb == nblks) i_endidx = npromz
      DO jc=1,i_endidx
        iidx = idx_old2new(idx_1d(jc,jb))
        jc_new = idx_no(iidx)
        jb_new = blk_no(iidx)
        arr(jc_new,jb_new) = tmp(jc,jb)
      END DO
    END DO
  END SUBROUTINE reorder_array_pos_g2D


  !> reorder array entries according to a given permutation array.
  !  2D implementation, t_tangent_vectors.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_t2D(arr, idx_old2new, nblks, npromz)
    TYPE(t_tangent_vectors), INTENT(INOUT), CONTIGUOUS :: arr(:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS :: idx_old2new(:) ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz

    ! local variables
    TYPE(t_tangent_vectors) :: tmp(SIZE(arr,1), SIZE(arr,2))
    INTEGER  :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:) = arr(:,:)
    DO jb=1,nblks
      i_endidx = nproma
      if (jb == nblks) i_endidx = npromz
      DO jc=1,i_endidx
        iidx = idx_old2new(idx_1d(jc,jb))
        jc_new = idx_no(iidx)
        jb_new = blk_no(iidx)
        arr(jc_new,jb_new) = tmp(jc,jb)
      END DO
    END DO
  END SUBROUTINE reorder_array_pos_t2D


  !> reorder array entries according to a given permutation array.
  !  2D implementation, t_cartesian_coordinates
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_c2D(arr, idx_old2new, nblks, npromz)
    TYPE(t_cartesian_coordinates), INTENT(INOUT), CONTIGUOUS :: arr(:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS :: idx_old2new(:) ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz

    ! local variables
    TYPE(t_cartesian_coordinates) :: tmp(SIZE(arr,1), SIZE(arr,2))
    INTEGER  :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:) = arr(:,:)
    DO jb=1,nblks
      i_endidx = nproma
      if (jb == nblks) i_endidx = npromz
      DO jc=1,i_endidx
        iidx = idx_old2new(idx_1d(jc,jb))
        jc_new = idx_no(iidx)
        jb_new = blk_no(iidx)
        arr(jc_new,jb_new) = tmp(jc,jb)
      END DO
    END DO
  END SUBROUTINE reorder_array_pos_c2D


  !> reorder array entries according to a given permutation array.
  !  3D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_i3D(arr, idx_old2new, nblks, npromz, &
    &                              idim_nproma, idim_blks, opt_owner)
    INTEGER,     INTENT(INOUT), CONTIGUOUS :: arr(:,:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz
    INTEGER,     INTENT(IN)    :: idim_nproma, idim_blks ! index positions of nproma and blocks
    LOGICAL, INTENT(IN), OPTIONAL, CONTIGUOUS :: opt_owner(:,:)

    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_array_pos_i3D"
    INTEGER :: tmp(SIZE(arr,1), SIZE(arr,2), SIZE(arr,3))
    INTEGER :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:,:) = arr(:,:,:)
    IF ((idim_nproma == 1) .AND. (idim_blks==2)) THEN
      IF (PRESENT(opt_owner)) THEN
        DO jb=1,nblks
          i_endidx = nproma
          IF (jb == nblks) i_endidx = npromz
          DO jc=1,i_endidx
            IF (opt_owner(jc,jb)) THEN
              iidx = idx_old2new(idx_1d(jc,jb))
              jc_new = idx_no(iidx)
              jb_new = blk_no(iidx)
              arr(jc_new,jb_new,:) = tmp(jc,jb,:)
            END IF
          END DO
        END DO
      ELSE
        DO jb=1,nblks
          i_endidx = nproma
          IF (jb == nblks) i_endidx = npromz
          DO jc=1,i_endidx
            iidx = idx_old2new(idx_1d(jc,jb))
            jc_new = idx_no(iidx)
            jb_new = blk_no(iidx)
            arr(jc_new,jb_new,:) = tmp(jc,jb,:)
          END DO
        END DO
      END IF
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
  END SUBROUTINE reorder_array_pos_i3D


  !> reorder array entries according to a given permutation array.
  !  4D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_i4D(arr, idx_old2new, nblks, npromz, &
    &                              idim_nproma, idim_blks)
    INTEGER,     INTENT(INOUT), CONTIGUOUS :: arr(:,:,:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz
    INTEGER,     INTENT(IN)    :: idim_nproma, idim_blks ! index positions of nproma and blocks

    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_array_pos_i4D"
    INTEGER :: tmp(SIZE(arr,1), SIZE(arr,2), SIZE(arr,3), SIZE(arr,4))
    INTEGER :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:,:,:) = arr(:,:,:,:)
    IF ((idim_nproma == 1) .AND. (idim_blks==2)) THEN
      DO jb=1,nblks
        i_endidx = nproma
        IF (jb == nblks) i_endidx = npromz
        DO jc=1,i_endidx
          iidx = idx_old2new(idx_1d(jc,jb))
          jc_new = idx_no(iidx)
          jb_new = blk_no(iidx)
          arr(jc_new,jb_new,:,:) = tmp(jc,jb,:,:)
        END DO
      END DO
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
  END SUBROUTINE reorder_array_pos_i4D


  !> reorder array entries according to a given permutation array.
  !  3D implementation, REALS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_r3D(arr, idx_old2new, nblks, npromz, &
    &                              idim_nproma, idim_blks)
    REAL(wp),    INTENT(INOUT), CONTIGUOUS :: arr(:,:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz
    INTEGER,     INTENT(IN)    :: idim_nproma, idim_blks ! index positions of nproma and blocks

    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_array_pos_r3D"
    REAL(wp) :: tmp(SIZE(arr,1), SIZE(arr,2), SIZE(arr,3))
    INTEGER  :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:,:) = arr(:,:,:)
    IF ((idim_nproma == 1) .AND. (idim_blks==2)) THEN
      DO jb=1,nblks
        i_endidx = nproma
        IF (jb == nblks) i_endidx = npromz
        DO jc=1,i_endidx
          iidx = idx_old2new(idx_1d(jc,jb))
          jc_new = idx_no(iidx)
          jb_new = blk_no(iidx)
          arr(jc_new,jb_new,:) = tmp(jc,jb,:)
        END DO
      END DO
    ELSE IF ((idim_nproma == 1) .AND. (idim_blks==3)) THEN
      DO jb=1,nblks
        i_endidx = nproma
        IF (jb == nblks) i_endidx = npromz
        DO jc=1,i_endidx
          iidx = idx_old2new(idx_1d(jc,jb))
          jc_new = idx_no(iidx)
          jb_new = blk_no(iidx)
          arr(jc_new,:,jb_new) = tmp(jc,:,jb)
        END DO
      END DO
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
  END SUBROUTINE reorder_array_pos_r3D


  !> reorder array entries according to a given permutation array.
  !  3D implementation, REALS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_pos_t3D(arr, idx_old2new, nblks, npromz, &
    &                              idim_nproma, idim_blks)
    TYPE(t_tangent_vectors),    INTENT(INOUT), CONTIGUOUS :: arr(:,:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz
    INTEGER,     INTENT(IN)    :: idim_nproma, idim_blks ! index positions of nproma and blocks

    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_array_pos_t3D"
    TYPE(t_tangent_vectors) :: tmp(SIZE(arr,1), SIZE(arr,2), SIZE(arr,3))
    INTEGER  :: iidx,jc,jb,jc_new,jb_new,i_endidx

    tmp(:,:,:) = arr(:,:,:)
    IF ((idim_nproma == 1) .AND. (idim_blks==2)) THEN
      DO jb=1,nblks
        i_endidx = nproma
        IF (jb == nblks) i_endidx = npromz
        DO jc=1,i_endidx
          iidx = idx_old2new(idx_1d(jc,jb))
          jc_new = idx_no(iidx)
          jb_new = blk_no(iidx)
          arr(jc_new,jb_new,:) = tmp(jc,jb,:)
        END DO
      END DO
    ELSE IF ((idim_nproma == 1) .AND. (idim_blks==3)) THEN
      DO jb=1,nblks
        i_endidx = nproma
        IF (jb == nblks) i_endidx = npromz
        DO jc=1,i_endidx
          iidx = idx_old2new(idx_1d(jc,jb))
          jc_new = idx_no(iidx)
          jb_new = blk_no(iidx)
          arr(jc_new,:,jb_new) = tmp(jc,:,jb)
        END DO
      END DO
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
  END SUBROUTINE reorder_array_pos_t3D


  !> translate array content according to a given permutation array.
  !  3D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_content_i1D(arr, glb_idx, idx_old2new, nentries)
    INTEGER,     INTENT(INOUT), CONTIGUOUS :: arr(:)                 ! global-size array
    INTEGER,     INTENT(IN), CONTIGUOUS    :: glb_idx(:)             ! index local->global
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nentries               ! no. of local entries

    ! local variables
    INTEGER :: j, idx_g
    INTEGER :: tmp(SIZE(arr,1))

    tmp(:) = arr(:)
    arr(:) = 0
    DO j=1,nentries
      idx_g = glb_idx(j)
      arr(idx_g) = idx_old2new(tmp(idx_g))
    END DO
  END SUBROUTINE reorder_array_content_i1D


  !> translate array content according to a given permutation array.
  !  3D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_content_i2D(idx_arr, blk_arr, idx_old2new, nblks, npromz, &
    &                                 opt_owner)
    INTEGER,     INTENT(INOUT), CONTIGUOUS :: idx_arr(:,:), blk_arr(:,:)
    INTEGER,     INTENT(IN), CONTIGUOUS    :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz
    LOGICAL, INTENT(IN), OPTIONAL, CONTIGUOUS :: opt_owner(:,:)

    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_array_content_i2D"
    INTEGER :: iidx,jc,jb,i_endidx

    DO jb=1,nblks
      i_endidx = MERGE(nproma, npromz, jb /= nblks)

      IF (PRESENT(opt_owner)) THEN
        DO jc=1,i_endidx
          IF (opt_owner(jc,jb)) THEN
            iidx = idx_old2new(idx_1d(idx_arr(jc,jb), blk_arr(jc,jb)))
            idx_arr(jc,jb) = idx_no(iidx)
            blk_arr(jc,jb) = blk_no(iidx)
          END IF
        END DO
      ELSE
        DO jc=1,i_endidx
          iidx = idx_old2new(idx_1d(idx_arr(jc,jb), blk_arr(jc,jb)))
          idx_arr(jc,jb) = idx_no(iidx)
          blk_arr(jc,jb) = blk_no(iidx)
        END DO
      END IF
    END DO
  END SUBROUTINE reorder_array_content_i2D


  !> translate array content according to a given permutation array.
  !  3D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_content_i3D(idx_arr, blk_arr, idx_old2new, nblks, npromz, &
    &                                  idim_nproma, idim_blks, opt_owner, opt_lcatch_zeros)
    INTEGER,     INTENT(INOUT) :: idx_arr(:,:,:), blk_arr(:,:,:)
    INTEGER,     INTENT(IN)    :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz
    INTEGER,     INTENT(IN)    :: idim_nproma, idim_blks ! index positions of nproma and blocks
    LOGICAL, INTENT(IN), OPTIONAL :: opt_owner(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_lcatch_zeros
    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_array_content_i3D"
    INTEGER :: iidx,jc,jb,jk,i_endidx,dim3
    LOGICAL :: lcatch_zeros

    lcatch_zeros = .TRUE.
    IF (PRESENT(opt_lcatch_zeros))  lcatch_zeros = opt_lcatch_zeros

    IF ((idim_nproma == 1) .AND. (idim_blks==2)) THEN
      dim3 = SIZE(idx_arr,3)
      DO jb=1,nblks
        i_endidx = nproma
        IF (jb == nblks) i_endidx = npromz

        IF (lcatch_zeros) THEN
          IF (PRESENT(opt_owner)) THEN
            DO jk=1,dim3
              DO jc=1,i_endidx
                IF (.NOT. opt_owner(jc,jb)) CYCLE
                IF (idx_arr(jc,jb,jk) <= 0) CYCLE
                iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk), blk_arr(jc,jb,jk)))
                idx_arr(jc,jb,jk) = idx_no(iidx)
                blk_arr(jc,jb,jk) = blk_no(iidx)
              END DO
            END DO
          ELSE
            DO jk=1,dim3
              DO jc=1,i_endidx
                IF (idx_arr(jc,jb,jk) <= 0) CYCLE
                iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk), blk_arr(jc,jb,jk)))
                idx_arr(jc,jb,jk) = idx_no(iidx)
                blk_arr(jc,jb,jk) = blk_no(iidx)
              END DO
            END DO
          END IF
        ELSE
          IF (PRESENT(opt_owner)) THEN
            DO jk=1,dim3
              DO jc=1,i_endidx
                IF (.NOT. opt_owner(jc,jb)) CYCLE
                iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk), blk_arr(jc,jb,jk)))
                idx_arr(jc,jb,jk) = idx_no(iidx)
                blk_arr(jc,jb,jk) = blk_no(iidx)
              END DO
            END DO
          ELSE
            DO jk=1,dim3
              DO jc=1,i_endidx
                iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk), blk_arr(jc,jb,jk)))
                idx_arr(jc,jb,jk) = idx_no(iidx)
                blk_arr(jc,jb,jk) = blk_no(iidx)
              END DO
            END DO
          END IF
        END IF
      END DO
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
  END SUBROUTINE reorder_array_content_i3D


  !> translate array content according to a given permutation array.
  !  3D implementation, INTEGERS.
  !
  !  @todo OpenMP parallelization.
  !
  SUBROUTINE reorder_array_content_i4D(idx_arr, blk_arr, idx_old2new, nblks, npromz, &
    &                                 idim_nproma, idim_blks, opt_owner, opt_lcatch_zeros)
    INTEGER,     INTENT(INOUT) :: idx_arr(:,:,:,:), blk_arr(:,:,:,:)
    INTEGER,     INTENT(IN)    :: idx_old2new(:)         ! permutation array
    INTEGER,     INTENT(IN)    :: nblks, npromz
    INTEGER,     INTENT(IN)    :: idim_nproma, idim_blks ! index positions of nproma and blocks
    LOGICAL, INTENT(IN), OPTIONAL :: opt_owner(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: opt_lcatch_zeros
    ! local variables
    CHARACTER (LEN=*), PARAMETER :: routine = modname//":reorder_array_content_i4D"
    INTEGER :: iidx,jc,jb,jk1,jk2,i_endidx,dim3,dim4
    LOGICAL :: lcatch_zeros

    lcatch_zeros = .FALSE.
    IF (PRESENT(opt_lcatch_zeros))  lcatch_zeros = opt_lcatch_zeros

    IF ((idim_nproma == 1) .AND. (idim_blks==2)) THEN
      dim3 = SIZE(idx_arr,3)
      dim4 = SIZE(idx_arr,4)
      DO jb=1,nblks
        i_endidx = nproma
        IF (jb == nblks) i_endidx = npromz

        IF (lcatch_zeros) THEN
          IF (PRESENT(opt_owner)) THEN
            DO jk1=1,dim3
              DO jk2=1,dim4
                DO jc=1,i_endidx
                  IF (.NOT. opt_owner(jc,jb)) CYCLE
                  IF (idx_arr(jc,jb,jk1,jk2) <= 0) CYCLE
                  iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk1,jk2), blk_arr(jc,jb,jk1,jk2)))
                  idx_arr(jc,jb,jk1,jk2) = idx_no(iidx)
                  blk_arr(jc,jb,jk1,jk2) = blk_no(iidx)
                END DO
              END DO
            END DO
          ELSE
            DO jk1=1,dim3
              DO jk2=1,dim4
                DO jc=1,i_endidx
                  IF (idx_arr(jc,jb,jk1,jk2) <= 0) CYCLE
                  iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk1,jk2), blk_arr(jc,jb,jk1,jk2)))
                  idx_arr(jc,jb,jk1,jk2) = idx_no(iidx)
                  blk_arr(jc,jb,jk1,jk2) = blk_no(iidx)
                END DO
              END DO
            END DO
          END IF
        ELSE
          IF (PRESENT(opt_owner)) THEN
            DO jk1=1,dim3
              DO jk2=1,dim4
                DO jc=1,i_endidx
                  IF (.NOT. opt_owner(jc,jb)) CYCLE
                  iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk1,jk2), blk_arr(jc,jb,jk1,jk2)))
                  idx_arr(jc,jb,jk1,jk2) = idx_no(iidx)
                  blk_arr(jc,jb,jk1,jk2) = blk_no(iidx)
                END DO
              END DO
            END DO
          ELSE
            DO jk1=1,dim3
              DO jk2=1,dim4
                DO jc=1,i_endidx
                  iidx = idx_old2new(idx_1d(idx_arr(jc,jb,jk1,jk2), blk_arr(jc,jb,jk1,jk2)))
                  idx_arr(jc,jb,jk1,jk2) = idx_no(iidx)
                  blk_arr(jc,jb,jk1,jk2) = blk_no(iidx)
                END DO
              END DO
            END DO
          END IF
        END IF
      END DO
    ELSE
      CALL finish(routine, "Internal error!")
    END IF
  END SUBROUTINE reorder_array_content_i4D

END MODULE mo_reorder_patches
