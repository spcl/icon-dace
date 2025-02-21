! The module <i>mo_alloc_domain</i>
! provides routines to allocate and deallocate patches. The routines
! included here have been extracted from mo_model_import_domain in order
! to resolve USE dependencies that led to large compile times on some
! platforms
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

MODULE mo_alloc_patches
  !-------------------------------------------------------------------------

  USE mo_kind,                    ONLY: wp
  USE mo_impl_constants,          ONLY: success, max_char_length,                    &
    &                                   min_rlcell, max_rlcell,                      &
    &                                   min_rledge, max_rledge,                      &
    &                                   min_rlvert, max_rlvert
  USE mo_exception,               ONLY: message, finish
  USE mo_model_domain,            ONLY: t_patch, t_pre_patch, c_num_edges,           &
    &                                   c_parent, c_child, c_phys_id, c_neighbor,    &
    &                                   c_edge, c_vertex, c_center,                  &
    &                                   c_refin_ctrl, e_parent, e_child, e_cell,     &
    &                                   e_refin_ctrl, v_cell, v_num_edges, v_vertex, &
    &                                   v_refin_ctrl
  USE mo_decomposition_tools,     ONLY: t_grid_domain_decomp_info,                   &
    &                                   init_glb2loc_index_lookup,                   &
    &                                   deallocate_glb2loc_index_lookup,             &
    &                                   uniform_partition,                           &
    &                                   partidx_of_elem_uniform_deco
  USE mo_parallel_config,         ONLY: nproma, num_dist_array_replicas
  USE mo_grid_config,             ONLY: n_dom_start, n_dom
  USE mo_mpi,                     ONLY: p_pe_work, p_comm_work, p_n_work



  USE mo_read_netcdf_distributed, ONLY: delete_distrib_read
  USE ppm_distributed_array,      ONLY: global_array_desc,                           &
    &                                   dist_mult_array_new,                         &
    &                                   dist_mult_array_delete,                      &
    &                                   ppm_int, ppm_real_dp



  USE ppm_extents,                ONLY: extent, extent_start, extent_size
  USE mo_communication,           ONLY: delete_comm_pattern,                         &
    &                                   delete_comm_gather_pattern

  IMPLICIT NONE

  PRIVATE

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: destruct_patches
  PUBLIC :: destruct_comm_patterns
  PUBLIC :: allocate_basic_patch
  PUBLIC :: allocate_pre_patch
  PUBLIC :: deallocate_pre_patch
  PUBLIC :: allocate_remaining_patch

  !-------------------------------------------------------------------------

CONTAINS


  !-------------------------------------------------------------------------
  !! Deallocation of all memory that was allocated
  !! to store the patch information of all levels.
  !!
  SUBROUTINE deallocate_patch( p_patch, lddmode )

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch
    LOGICAL,     OPTIONAL, INTENT(in)    :: lddmode  ! Skip fields not allocated in ldump_dd mode

    !local variables
    CHARACTER(*), PARAMETER :: routine ="mo_alloc_patches:deallocate_patch"
    INTEGER :: ist
    LOGICAL :: l_ddmode
    !-----------------------------------------------------------------------
    !
    IF (PRESENT(lddmode)) THEN
      l_ddmode = lddmode
    ELSE
      l_ddmode = .FALSE.
    ENDIF

    ! Deallocate grid information
    !
    ! Include deallocation of basic patch variables
! GZ: this directive should not be needed here, but the NEC compiler crashes otherwise (internal error in optimization phase)
!CDIR NOIEXPAND
    CALL deallocate_basic_patch(p_patch)

    call delete_distrib_read(p_patch%cells%dist_io_data)
    call delete_distrib_read(p_patch%edges%dist_io_data)
    call delete_distrib_read(p_patch%verts%dist_io_data)

    DEALLOCATE(p_patch%cells%dist_io_data, p_patch%edges%dist_io_data, &
      &        p_patch%verts%dist_io_data)

    call deallocate_decomp_info(p_patch%cells%decomp_info)
    call deallocate_decomp_info(p_patch%edges%decomp_info)
    call deallocate_decomp_info(p_patch%verts%decomp_info)

    DEALLOCATE(p_patch%cells%phys_id,    stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell physical ID failed')
    ENDIF

    IF (l_ddmode) RETURN

    DEALLOCATE( p_patch%cells%edge_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell edge_orientation failed')
    ENDIF
    DEALLOCATE( p_patch%cells%area,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell area failed')
    ENDIF
    DEALLOCATE( p_patch%cells%f_c,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch cell f_c failed')
    ENDIF
    !
    !
    DEALLOCATE(p_patch%edges%phys_id,    stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge physical ID failed')
    ENDIF
    DEALLOCATE( p_patch%edges%cell_idx,  &
      & p_patch%edges%cell_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge cell index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%vertex_idx,  &
      & p_patch%edges%vertex_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge vertex index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%tangent_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge system orientation failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_idx,  &
      & p_patch%edges%quad_blk,  stat=ist)
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge quad index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge quad orientation failed')
    ENDIF
    DEALLOCATE( p_patch%edges%butterfly_idx,  &
      & p_patch%edges%butterfly_blk,  stat=ist)
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge butterfly index failed')
    ENDIF
    DEALLOCATE( p_patch%edges%center,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge center failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal_normal failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual_normal failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal_cell,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal_normal_cell failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal_cell,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual_normal_cell failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_normal_vert,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal_normal_vert failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_normal_vert,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual_normal_vert failed')
    ENDIF
    DEALLOCATE( p_patch%edges%primal_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge primal edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_primal_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for inverse patch edge primal edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%dual_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge dual edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_dual_edge_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for inverse patch edge dual edge length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%edge_vert_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for edge_vert_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%inv_vert_vert_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for inverse vert_vert_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%edge_cell_length,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for edge_cell_length failed')
    ENDIF
    DEALLOCATE( p_patch%edges%area_edge,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge area_edge failed')
    ENDIF
    DEALLOCATE( p_patch%edges%quad_area,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge quad area failed')
    ENDIF
    DEALLOCATE( p_patch%edges%f_e,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge f_e failed')
    ENDIF
    DEALLOCATE( p_patch%edges%fn_e,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge fn_e failed')
    ENDIF
    DEALLOCATE( p_patch%edges%ft_e,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch edge ft_e failed')
    ENDIF
    !
    !
    DEALLOCATE( p_patch%verts%neighbor_idx,  &
      & p_patch%verts%neighbor_blk,  &
      & p_patch%verts%phys_id,       stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex neighbor index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%cell_idx,  &
      & p_patch%verts%cell_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex cell index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%edge_idx,  &
      & p_patch%verts%edge_blk,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex edge index failed')
    ENDIF
    DEALLOCATE( p_patch%verts%edge_orientation,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex edge orientation failed')
    ENDIF
    DEALLOCATE( p_patch%verts%num_edges,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex number of edges failed')
    ENDIF
    DEALLOCATE( p_patch%verts%dual_area,  stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex dual area failed')
    ENDIF
    DEALLOCATE( p_patch%verts%f_v, stat=ist )
    IF(ist/=success)THEN
      CALL finish  (routine,  'deallocate for patch vertex f_v failed')
    ENDIF


    CALL deallocate_patch_cartesian( p_patch )

  CONTAINS

    SUBROUTINE deallocate_decomp_info( decomp_info )

      TYPE(t_grid_domain_decomp_info), INTENT(inout) :: decomp_info

      INTEGER :: ist

      DEALLOCATE(decomp_info%glb_index, &
        &        decomp_info%owner_local, &
        &        decomp_info%owner_mask, &
        &        decomp_info%decomp_domain, stat=ist )
      IF(ist/=success)THEN
        CALL finish  (routine,  'deallocate in deallocate_decomp_info failed')
      ENDIF
      CALL deallocate_glb2loc_index_lookup(decomp_info%glb2loc_index)
    END SUBROUTINE deallocate_decomp_info

  END SUBROUTINE deallocate_patch

  !-------------------------------------------------------------------------

  SUBROUTINE allocate_patch_cartesian( p_patch )
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch

    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_alloc_patches:allocate_patch_cartesian'

    ALLOCATE( p_patch%edges%primal_cart_normal(nproma,p_patch%nblks_e) )
    p_patch%edges%primal_cart_normal(:,:)%x(1) = 0._wp
    p_patch%edges%primal_cart_normal(:,:)%x(2) = 0._wp
    p_patch%edges%primal_cart_normal(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%edges%dual_cart_normal(nproma,p_patch%nblks_e) )
    p_patch%edges%dual_cart_normal(:,:)%x(1) = 0._wp
    p_patch%edges%dual_cart_normal(:,:)%x(2) = 0._wp
    p_patch%edges%dual_cart_normal(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%cells%cartesian_center(nproma,p_patch%alloc_cell_blocks) )
    p_patch%cells%cartesian_center(:,:)%x(1) = 0._wp
    p_patch%cells%cartesian_center(:,:)%x(2) = 0._wp
    p_patch%cells%cartesian_center(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%edges%cartesian_center(nproma,p_patch%nblks_e) )
    p_patch%edges%cartesian_center(:,:)%x(1) = 0._wp
    p_patch%edges%cartesian_center(:,:)%x(2) = 0._wp
    p_patch%edges%cartesian_center(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%edges%cartesian_dual_middle(nproma,p_patch%nblks_e) )
    p_patch%edges%cartesian_dual_middle(:,:)%x(1) = 0._wp
    p_patch%edges%cartesian_dual_middle(:,:)%x(2) = 0._wp
    p_patch%edges%cartesian_dual_middle(:,:)%x(3) = 0._wp

    ALLOCATE( p_patch%verts%cartesian(nproma,p_patch%nblks_v) )
    p_patch%verts%cartesian(:,:)%x(1) = 0._wp
    p_patch%verts%cartesian(:,:)%x(2) = 0._wp
    p_patch%verts%cartesian(:,:)%x(3) = 0._wp

  END SUBROUTINE allocate_patch_cartesian
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  SUBROUTINE deallocate_patch_cartesian( p_patch )
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch

    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: &
      & method_name = 'mo_alloc_patches:deallocate_patch_cartesian'

    DEALLOCATE( p_patch%edges%primal_cart_normal, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%primal_cart_normal)")
    ENDIF

    DEALLOCATE( p_patch%edges%dual_cart_normal, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%dual_cart_normal)")
    ENDIF

    DEALLOCATE( p_patch%cells%cartesian_center, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%cells%cartesian_center)")
    ENDIF

    DEALLOCATE( p_patch%edges%cartesian_center, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%cartesian_center)")
    ENDIF

    DEALLOCATE( p_patch%edges%cartesian_dual_middle, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%edges%cartesian_dual_middle)")
    ENDIF

    DEALLOCATE( p_patch%verts%cartesian, stat=return_status )
    IF ( return_status /= success) THEN
      CALL finish  (method_name,"DEALLOCATE( p_patch%verts%cartesian)")
    ENDIF

  END SUBROUTINE deallocate_patch_cartesian
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> destruct_patches: Calls destruct_patch in a loop
  SUBROUTINE destruct_patches( p_patch )

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch(n_dom_start:)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_alloc_patches:destruct_patches'

    !local variables
    INTEGER :: jg, istart, iend
    !-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    istart = LBOUND(p_patch, 1)
    iend   = UBOUND(p_patch, 1)

    grid_level_loop: DO jg = istart, iend
      CALL deallocate_patch(p_patch(jg))
    END DO grid_level_loop

    CALL message (routine, 'destruct_patches finished')

  END SUBROUTINE destruct_patches
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> destruct_comm_patterns
  SUBROUTINE destruct_comm_patterns( p_patch, p_patch_local_parent )

    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch(n_dom_start:)
    TYPE(t_patch), TARGET, INTENT(inout) :: p_patch_local_parent(n_dom_start+1:)

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_alloc_patches:destruct_comm_patterns'

    !local variables
    INTEGER :: jg
    !-----------------------------------------------------------------------

    CALL message (TRIM(routine), 'start')

    DO jg = n_dom_start, n_dom

      CALL delete_comm_pattern(p_patch(jg)%comm_pat_c)
      CALL delete_comm_pattern(p_patch(jg)%comm_pat_c1)
      CALL delete_comm_pattern(p_patch(jg)%comm_pat_v)
      CALL delete_comm_pattern(p_patch(jg)%comm_pat_e)

      CALL delete_comm_gather_pattern(p_patch(jg)%comm_pat_gather_c)
      CALL delete_comm_gather_pattern(p_patch(jg)%comm_pat_gather_v)
      CALL delete_comm_gather_pattern(p_patch(jg)%comm_pat_gather_e)

      IF (jg > n_dom_start) THEN

        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_c)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_v)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_e)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_c1)

        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_glb_to_loc_c)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_glb_to_loc_e)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_loc_to_glb_c_fbk)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_loc_to_glb_e_fbk)
      ENDIF

    ENDDO

    CALL message (routine, 'destruct_comm_patterns finished')

  END SUBROUTINE destruct_comm_patterns
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !! Allocates all arrays in a basic patch
  !! These are the arrays which have to be read in full size on every PE
  !! since they are needed for domain decomposition.
  !!
  SUBROUTINE allocate_basic_patch(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch

    INTEGER :: max_childdom, max_cell_connectivity, max_vertex_connectivity

    ! Please note: The following variables in the patch MUST already be set:
    ! - alloc_cell_blocks
    ! - nblks_e
    ! - nblks_v
    ! - n_patch_cells
    ! - n_patch_edges
    ! - n_patch_verts
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! - max_childdom

    max_cell_connectivity = p_patch%cells%max_connectivity
    max_vertex_connectivity = p_patch%verts%max_connectivity
    p_patch%geometry_info%cell_type = p_patch%cells%max_connectivity ! this should be read by the grid !
    max_childdom = p_patch%max_childdom
    !
    ! !grid cells
    !
    p_patch%cells%dummy_cell_block = 0
    p_patch%cells%dummy_cell_index = 0
    ALLOCATE( p_patch%cells%num_edges(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%parent_loc_idx(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%parent_loc_blk(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%parent_glb_idx(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%parent_glb_blk(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%pc_idx(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%child_idx(nproma,p_patch%alloc_cell_blocks,4) )
    ALLOCATE( p_patch%cells%child_blk(nproma,p_patch%alloc_cell_blocks,4) )
    ALLOCATE( p_patch%cells%child_id(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%phys_id(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%neighbor_idx(nproma,p_patch%alloc_cell_blocks,max_cell_connectivity) )
    ALLOCATE( p_patch%cells%neighbor_blk(nproma,p_patch%alloc_cell_blocks,max_cell_connectivity) )
    ALLOCATE( p_patch%cells%edge_idx(nproma,p_patch%alloc_cell_blocks,max_cell_connectivity) )
    ALLOCATE( p_patch%cells%edge_blk(nproma,p_patch%alloc_cell_blocks,max_cell_connectivity) )
    ALLOCATE( p_patch%cells%vertex_idx(nproma,p_patch%alloc_cell_blocks,max_cell_connectivity) )
    ALLOCATE( p_patch%cells%vertex_blk(nproma,p_patch%alloc_cell_blocks,max_cell_connectivity) )
    ALLOCATE( p_patch%cells%center(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%refin_ctrl(nproma,p_patch%alloc_cell_blocks) )
    ALLOCATE( p_patch%cells%start_idx(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%end_idx(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%start_blk(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%end_blk(min_rlcell:max_rlcell,max_childdom) )
    ALLOCATE( p_patch%cells%start_index(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch%cells%end_index(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch%cells%start_block(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch%cells%end_block(min_rlcell:max_rlcell) )

    !
    ! !grid edges
    !
    ALLOCATE( p_patch%edges%parent_loc_idx(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%parent_loc_blk(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%parent_glb_idx(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%parent_glb_blk(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%pc_idx(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%child_idx(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%child_blk(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%child_id(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%refin_ctrl(nproma,p_patch%nblks_e) )
    ALLOCATE( p_patch%edges%start_idx(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%end_idx(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%start_blk(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%end_blk(min_rledge:max_rledge,max_childdom) )
    ALLOCATE( p_patch%edges%cell_idx(nproma,p_patch%nblks_e,2) )
    ALLOCATE( p_patch%edges%cell_blk(nproma,p_patch%nblks_e,2) )
    ALLOCATE( p_patch%edges%vertex_idx(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%vertex_blk(nproma,p_patch%nblks_e,4) )
    ALLOCATE( p_patch%edges%start_index(min_rledge:max_rledge) )
    ALLOCATE( p_patch%edges%end_index(min_rledge:max_rledge) )
    ALLOCATE( p_patch%edges%start_block(min_rledge:max_rledge) )
    ALLOCATE( p_patch%edges%end_block(min_rledge:max_rledge) )

    !
    ! !grid verts
    !
    ALLOCATE( p_patch%verts%vertex(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%refin_ctrl(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%start_idx(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%end_idx(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%start_blk(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%end_blk(min_rlvert:max_rlvert,max_childdom) )
    ALLOCATE( p_patch%verts%cell_idx(nproma,p_patch%nblks_v,max_vertex_connectivity) )
    ALLOCATE( p_patch%verts%cell_blk(nproma,p_patch%nblks_v,max_vertex_connectivity) )
    ALLOCATE( p_patch%verts%num_edges(nproma,p_patch%nblks_v) )
    ALLOCATE( p_patch%verts%start_index(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch%verts%end_index(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch%verts%start_block(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch%verts%end_block(min_rlvert:max_rlvert) )
    ! Set all newly allocated arrays to 0

    p_patch%cells%num_edges = 0
    p_patch%cells%parent_loc_idx = 0
    p_patch%cells%parent_loc_blk = 0
    p_patch%cells%parent_glb_idx = 0
    p_patch%cells%parent_glb_blk = 0
    p_patch%cells%pc_idx = 0
    p_patch%cells%child_idx = 0
    p_patch%cells%child_blk = 0
    p_patch%cells%child_id = 0
    p_patch%cells%phys_id = 0
    p_patch%cells%neighbor_idx = 0
    p_patch%cells%neighbor_blk = 0
    p_patch%cells%edge_idx = 0
    p_patch%cells%edge_blk = 0
    p_patch%cells%vertex_idx = 0
    p_patch%cells%vertex_blk = 0
    p_patch%cells%center(:,:)%lon = 0._wp
    p_patch%cells%center(:,:)%lat = 0._wp
    p_patch%cells%refin_ctrl = 0
    p_patch%cells%start_idx = 0
    p_patch%cells%end_idx = 0
    p_patch%cells%start_blk = 0
    p_patch%cells%end_blk = 0
    p_patch%cells%start_index = 0
    p_patch%cells%end_index = 0
    p_patch%cells%start_block = 0
    p_patch%cells%end_block = 0

    p_patch%edges%parent_loc_idx = 0
    p_patch%edges%parent_loc_blk = 0
    p_patch%edges%parent_glb_idx = 0
    p_patch%edges%parent_glb_blk = 0
    p_patch%edges%pc_idx = 0
    p_patch%edges%child_idx = 0
    p_patch%edges%child_blk = 0
    p_patch%edges%child_id = 0
    p_patch%edges%refin_ctrl = 0
    p_patch%edges%start_idx = 0
    p_patch%edges%end_idx = 0
    p_patch%edges%start_blk = 0
    p_patch%edges%end_blk = 0
    p_patch%edges%start_index = 0
    p_patch%edges%end_index = 0
    p_patch%edges%start_block = 0
    p_patch%edges%end_block = 0

    p_patch%verts%vertex(:,:)%lon = 0._wp
    p_patch%verts%vertex(:,:)%lat = 0._wp
    p_patch%verts%refin_ctrl = 0
    p_patch%verts%start_idx = 0
    p_patch%verts%end_idx = 0
    p_patch%verts%start_blk = 0
    p_patch%verts%end_blk = 0
    p_patch%verts%cell_idx = 0
    p_patch%verts%cell_blk = 0
    p_patch%verts%num_edges = 0
    p_patch%verts%start_index = 0
    p_patch%verts%end_index = 0
    p_patch%verts%start_block = 0
    p_patch%verts%end_block = 0

  END SUBROUTINE allocate_basic_patch

  !-------------------------------------------------------------------------
  !! Allocates all arrays in a basic patch
  !! These are the arrays which have to be read in full size on every PE
  !! since they are needed for domain decomposition.
  !!
  SUBROUTINE allocate_pre_patch(p_patch_pre)

    TYPE(t_pre_patch), INTENT(inout) :: p_patch_pre

    INTEGER :: max_childdom
    TYPE(global_array_desc) :: dist_cell_desc(9), dist_edge_desc(4), &
         dist_vert_desc(4)

    TYPE(extent) :: &
      &             local_cell_chunks(2, 9), &
      &             local_edge_chunks(2, 4), &
      &             local_vert_chunks(2, 4), &
      &             process_space
    INTEGER :: num_replicas, replica_idx, mpierr, dist_array_comm, &
      &        dist_array_pes_start, dist_array_pes_size
    INTEGER :: max_cell_connectivity, max_vertex_connectivity

    ! Please note: The following variables in the patch MUST already be set:
    ! - alloc_cell_blocks
    ! - nblks_e
    ! - nblks_v
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! - max_childdom
    ! - cells%max_connectivity
    ! - verts%max_connectivity
    max_cell_connectivity = p_patch_pre%cells%max_connectivity
    max_vertex_connectivity = p_patch_pre%verts%max_connectivity
    p_patch_pre%cell_type = p_patch_pre%cells%max_connectivity ! this should be read by the grid !
    max_childdom = p_patch_pre%max_childdom

    ! some preliminary computation for the distributed data

    dist_cell_desc(c_num_edges)%a_rank = 1
    dist_cell_desc(c_num_edges)%rect(1)%first = 1
    dist_cell_desc(c_num_edges)%rect(1)%size = p_patch_pre%n_patch_cells_g
    dist_cell_desc(c_num_edges)%element_dt = ppm_int

    dist_cell_desc(c_parent)%a_rank = 1
    dist_cell_desc(c_parent)%rect(1)%first = 1
    dist_cell_desc(c_parent)%rect(1)%size = p_patch_pre%n_patch_cells_g
    dist_cell_desc(c_parent)%element_dt = ppm_int

    dist_cell_desc(c_child) = dist_cell_desc(c_num_edges)
    dist_cell_desc(c_child)%a_rank = 2
    dist_cell_desc(c_child)%rect(2) = extent(first=1, size = 4)

    dist_cell_desc(c_phys_id) = dist_cell_desc(c_num_edges)

    dist_cell_desc(c_neighbor) = dist_cell_desc(c_num_edges)
    dist_cell_desc(c_neighbor)%a_rank = 2
    dist_cell_desc(c_neighbor)%rect(2) &
         = extent(first=1, size = max_cell_connectivity)

    dist_cell_desc(c_edge) = dist_cell_desc(c_neighbor)

    dist_cell_desc(c_vertex) = dist_cell_desc(c_neighbor)

    dist_cell_desc(c_center)%a_rank = 2
    dist_cell_desc(c_center)%rect(1) &
         = extent(first = 1, size = p_patch_pre%n_patch_cells_g)
    dist_cell_desc(c_center)%rect(2) &
         = extent(first = 1, size = 2)
    dist_cell_desc(c_center)%element_dt = ppm_real_dp

    dist_cell_desc(c_refin_ctrl) = dist_cell_desc(c_num_edges)

    dist_edge_desc(e_parent)%a_rank = 1
    dist_edge_desc(e_parent)%rect(1)%first = 1
    dist_edge_desc(e_parent)%rect(1)%size = p_patch_pre%n_patch_edges_g
    dist_edge_desc(e_parent)%element_dt = ppm_int

    dist_edge_desc(e_child) = dist_edge_desc(e_parent)
    dist_edge_desc(e_child)%a_rank = 2
    dist_edge_desc(e_child)%rect(2) = extent(first = 1, size = 4)

    dist_edge_desc(e_cell) = dist_edge_desc(e_parent)
    dist_edge_desc(e_cell)%a_rank = 2
    dist_edge_desc(e_cell)%rect(2)%first = 1
    dist_edge_desc(e_cell)%rect(2)%size = 2

    dist_edge_desc(e_refin_ctrl) = dist_edge_desc(e_parent)

    dist_vert_desc(v_cell)%a_rank = 2
    dist_vert_desc(v_cell)%rect(1) &
         = extent(first = 1, size = p_patch_pre%n_patch_verts_g)
    dist_vert_desc(v_cell)%rect(2) = extent(first = 1, size = max_vertex_connectivity)
    dist_vert_desc(v_cell)%element_dt = ppm_int

    dist_vert_desc(v_num_edges)%a_rank = 1
    dist_vert_desc(v_num_edges)%rect(1)%first = 1
    dist_vert_desc(v_num_edges)%rect(1)%size = p_patch_pre%n_patch_verts_g
    dist_vert_desc(v_num_edges)%element_dt = ppm_int

    dist_vert_desc(v_vertex)%a_rank = 2
    dist_vert_desc(v_vertex)%rect(1) &
         = extent(first = 1, size = p_patch_pre%n_patch_verts_g)
    dist_vert_desc(v_vertex)%rect(2) &
         = extent(first = 1, size = 2)
    dist_vert_desc(v_vertex)%element_dt = ppm_real_dp

    dist_vert_desc(v_refin_ctrl) = dist_vert_desc(v_num_edges)

    process_space = extent(1, p_n_work)
    num_replicas = MAX(1, MIN(num_dist_array_replicas, p_n_work))
    replica_idx = partidx_of_elem_uniform_deco(process_space, num_replicas, &
      &                                        p_pe_work+1)




    dist_array_comm = p_comm_work

    p_patch_pre%dist_array_pes = &
      uniform_partition(process_space, num_replicas, replica_idx)
    p_patch_pre%dist_array_comm = dist_array_comm

    dist_array_pes_start = extent_start(p_patch_pre%dist_array_pes) - 1
    dist_array_pes_size = extent_size(p_patch_pre%dist_array_pes)

    p_patch_pre%cells%local_chunk(1,1) = &
      uniform_partition(dist_cell_desc(c_num_edges)%rect(1), &
        &               dist_array_pes_size, &
        &               p_pe_work + 1 - dist_array_pes_start)
    p_patch_pre%edges%local_chunk(1,1) = &
      uniform_partition(dist_edge_desc(e_parent)%rect(1), &
        &               dist_array_pes_size, &
        &               p_pe_work + 1 - dist_array_pes_start)
    p_patch_pre%verts%local_chunk(1,1) = &
      uniform_partition(dist_vert_desc(v_cell)%rect(1), &
        &               dist_array_pes_size, &
        &               p_pe_work + 1 - dist_array_pes_start)

    !
    ! !grid cells
    !
    local_cell_chunks(1, :) = p_patch_pre%cells%local_chunk(1, 1)
    local_cell_chunks(2, c_child) = extent(first=1, size=4)
    local_cell_chunks(2, c_neighbor:c_vertex) &
         = extent(first=1, size=max_cell_connectivity)
    local_cell_chunks(1, c_center) = p_patch_pre%cells%local_chunk(1,1)
    local_cell_chunks(2, c_center) = extent(first = 1, size = 2)

    p_patch_pre%cells%dist = dist_mult_array_new( &
         dist_cell_desc, local_cell_chunks, dist_array_comm, &



         cache_size=MIN(10, p_n_work) &

         )
    ALLOCATE( p_patch_pre%cells%start(min_rlcell:max_rlcell) )
    ALLOCATE( p_patch_pre%cells%end(min_rlcell:max_rlcell) )

    !
    ! !grid edges
    !
    local_edge_chunks(1, e_parent) = p_patch_pre%edges%local_chunk(1, 1)
    local_edge_chunks(1, e_child) = p_patch_pre%edges%local_chunk(1, 1)
    local_edge_chunks(2, e_child) = extent(first = 1, size = 4)
    local_edge_chunks(1, e_cell) = p_patch_pre%edges%local_chunk(1, 1)
    local_edge_chunks(2, e_cell) = extent(first = 1, size = 2)
    local_edge_chunks(1, e_refin_ctrl) = p_patch_pre%edges%local_chunk(1, 1)

    p_patch_pre%edges%dist = dist_mult_array_new( &
      dist_edge_desc, local_edge_chunks, dist_array_comm, &



         cache_size=MIN(10, p_n_work) &

         )
    ALLOCATE( p_patch_pre%edges%start(min_rledge:max_rledge) )
    ALLOCATE( p_patch_pre%edges%end(min_rledge:max_rledge) )

    !
    ! !grid verts
    !
    local_vert_chunks(1, v_cell) = p_patch_pre%verts%local_chunk(1,1)
    local_vert_chunks(2, v_cell) = extent(first = 1, size = max_vertex_connectivity)
    local_vert_chunks(1, v_num_edges) = p_patch_pre%verts%local_chunk(1,1)
    local_vert_chunks(1, v_vertex) = p_patch_pre%verts%local_chunk(1,1)
    local_vert_chunks(2, v_vertex) = extent(first = 1, size = 2)
    local_vert_chunks(1, v_refin_ctrl) = p_patch_pre%verts%local_chunk(1,1)

    p_patch_pre%verts%dist = dist_mult_array_new( &
      dist_vert_desc, local_vert_chunks, dist_array_comm, &



         cache_size=MIN(10, p_n_work) &

         )
    ALLOCATE( p_patch_pre%verts%start(min_rlvert:max_rlvert) )
    ALLOCATE( p_patch_pre%verts%end(min_rlvert:max_rlvert) )
    ! Set all newly allocated arrays to 0

    p_patch_pre%cells%start = 0
    p_patch_pre%cells%end = 0

    p_patch_pre%edges%start = 0
    p_patch_pre%edges%end = 0

    p_patch_pre%verts%start = 0
    p_patch_pre%verts%end = 0

  END SUBROUTINE allocate_pre_patch
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Deallocates all arrays in a basic patch
  !!
  SUBROUTINE deallocate_basic_patch(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch
    !
    ! !grid cells
    !
    DEALLOCATE( p_patch%cells%num_edges )
    IF (ALLOCATED(p_patch%cells%parent_loc_idx)) &
      DEALLOCATE( p_patch%cells%parent_loc_idx )
    IF (ALLOCATED(p_patch%cells%parent_loc_blk)) &
      DEALLOCATE( p_patch%cells%parent_loc_blk )
    IF (ALLOCATED(p_patch%cells%parent_glb_idx)) &
      DEALLOCATE( p_patch%cells%parent_glb_idx )
    IF (ALLOCATED(p_patch%cells%parent_glb_blk)) &
      DEALLOCATE( p_patch%cells%parent_glb_blk )
    DEALLOCATE( p_patch%cells%pc_idx )
    IF (ALLOCATED(p_patch%cells%child_idx)) &
      DEALLOCATE( p_patch%cells%child_idx )
    IF (ALLOCATED(p_patch%cells%child_blk)) &
      DEALLOCATE( p_patch%cells%child_blk )
    DEALLOCATE( p_patch%cells%child_id )
    DEALLOCATE( p_patch%cells%neighbor_idx )
    DEALLOCATE( p_patch%cells%neighbor_blk )
    DEALLOCATE( p_patch%cells%edge_idx )
    DEALLOCATE( p_patch%cells%edge_blk )
    DEALLOCATE( p_patch%cells%vertex_idx )
    DEALLOCATE( p_patch%cells%vertex_blk )
    DEALLOCATE( p_patch%cells%center )
    DEALLOCATE( p_patch%cells%refin_ctrl )
    DEALLOCATE( p_patch%cells%start_idx )
    DEALLOCATE( p_patch%cells%end_idx )
    DEALLOCATE( p_patch%cells%start_blk )
    DEALLOCATE( p_patch%cells%end_blk )
    DEALLOCATE( p_patch%cells%start_index )
    DEALLOCATE( p_patch%cells%end_index )
    DEALLOCATE( p_patch%cells%start_block )
    DEALLOCATE( p_patch%cells%end_block )
    !
    ! !grid edges
    !
    IF (ALLOCATED(p_patch%edges%parent_loc_idx)) &
      DEALLOCATE( p_patch%edges%parent_loc_idx )
    IF (ALLOCATED(p_patch%edges%parent_loc_blk)) &
      DEALLOCATE( p_patch%edges%parent_loc_blk )
    IF (ALLOCATED(p_patch%edges%parent_glb_idx)) &
      DEALLOCATE( p_patch%edges%parent_glb_idx )
    IF (ALLOCATED(p_patch%edges%parent_glb_blk)) &
      DEALLOCATE( p_patch%edges%parent_glb_blk )
    DEALLOCATE( p_patch%edges%pc_idx )
    IF (ALLOCATED(p_patch%edges%child_idx)) &
      DEALLOCATE( p_patch%edges%child_idx )
    IF (ALLOCATED(p_patch%edges%child_blk)) &
      DEALLOCATE( p_patch%edges%child_blk )
    DEALLOCATE( p_patch%edges%child_id )
    DEALLOCATE( p_patch%edges%refin_ctrl )
    DEALLOCATE( p_patch%edges%start_idx )
    DEALLOCATE( p_patch%edges%end_idx )
    DEALLOCATE( p_patch%edges%start_blk )
    DEALLOCATE( p_patch%edges%end_blk )
    DEALLOCATE( p_patch%edges%start_index )
    DEALLOCATE( p_patch%edges%end_index )
    DEALLOCATE( p_patch%edges%start_block )
    DEALLOCATE( p_patch%edges%end_block )
    ! DEALLOCATE( p_patch%edges%cell_idx )
    ! DEALLOCATE( p_patch%edges%cell_blk )
    ! DEALLOCATE( p_patch%edges%vertex_idx )
    ! DEALLOCATE( p_patch%edges%vertex_blk )
    !
    ! !grid verts
    !
    DEALLOCATE( p_patch%verts%vertex )
    DEALLOCATE( p_patch%verts%refin_ctrl )
    DEALLOCATE( p_patch%verts%start_idx )
    DEALLOCATE( p_patch%verts%end_idx )
    DEALLOCATE( p_patch%verts%start_blk )
    DEALLOCATE( p_patch%verts%end_blk )
    DEALLOCATE( p_patch%verts%start_index )
    DEALLOCATE( p_patch%verts%end_index )
    DEALLOCATE( p_patch%verts%start_block )
    DEALLOCATE( p_patch%verts%end_block )
    ! DEALLOCATE( p_patch%verts%cell_idx )
    ! DEALLOCATE( p_patch%verts%cell_blk )
    ! DEALLOCATE( p_patch%verts%num_edges )

  END SUBROUTINE deallocate_basic_patch

  !-------------------------------------------------------------------------
  !! Deallocates all arrays in a basic patch
  !!
  SUBROUTINE deallocate_pre_patch(p_patch_pre)

    TYPE(t_pre_patch), INTENT(inout) :: p_patch_pre
    INTEGER :: mpierr
    !
    ! !grid cells
    !
    CALL dist_mult_array_delete(p_patch_pre%cells%dist)
    DEALLOCATE( p_patch_pre%cells%start )
    DEALLOCATE( p_patch_pre%cells%end )
    !
    ! !grid edges
    !
    CALL dist_mult_array_delete(p_patch_pre%edges%dist)
    DEALLOCATE( p_patch_pre%edges%start )
    DEALLOCATE( p_patch_pre%edges%end )
    !
    ! !grid verts
    !
    CALL dist_mult_array_delete(p_patch_pre%verts%dist)
    DEALLOCATE( p_patch_pre%verts%start )
    DEALLOCATE( p_patch_pre%verts%end )



  END SUBROUTINE deallocate_pre_patch
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Allocates remaining arrays in a full patch.
  !! These are the arrays which have not been allocated in allocate_basic_patch
  !!
  SUBROUTINE allocate_remaining_patch(p_patch,iopmode)

    TYPE(t_patch), INTENT(inout) :: p_patch
    INTEGER,       INTENT(in)    :: iopmode

    INTEGER :: max_cell_connectivity, max_vertex_connectivity
    ! operation modes: iopmode=1: allocate all remaining patch arrays
    !                  iopmode=2: allocate only parallelization-related arrays
    !                  iopmode=3: allocate all but the parallelization-related arrays

    ! Please note: The following variables in the patch MUST already be set:
    ! - alloc_cell_blocks
    ! - nblks_e
    ! - nblks_v
    ! - n_patch_cells
    ! - n_patch_edges
    ! - n_patch_verts
    ! - n_patch_cells_g
    ! - n_patch_edges_g
    ! - n_patch_verts_g
    ! - max_childdom
    ! - cells%max_connectivity
    ! - verts%max_connectivity
    max_cell_connectivity   = p_patch%cells%max_connectivity
    max_vertex_connectivity = p_patch%verts%max_connectivity
    !
    ! !grid cells
    !
    IF (iopmode /= 2) THEN
      IF (.NOT. ALLOCATED(p_patch%cells%phys_id)) ALLOCATE( p_patch%cells%phys_id(nproma,p_patch%alloc_cell_blocks) )
      ALLOCATE( p_patch%cells%edge_orientation(nproma,p_patch%alloc_cell_blocks,max_cell_connectivity) )
      ALLOCATE( p_patch%cells%area(nproma,p_patch%alloc_cell_blocks) )
      ALLOCATE( p_patch%cells%f_c(nproma,p_patch%alloc_cell_blocks) )
    ENDIF

    !
    ! !grid edges
    !
    IF (iopmode /= 2) THEN
      ALLOCATE( p_patch%edges%phys_id(nproma,p_patch%nblks_e) )
      IF (ALLOCATED(p_patch%edges%cell_idx)) DEALLOCATE(p_patch%edges%cell_idx)
      ALLOCATE( p_patch%edges%cell_idx(nproma,p_patch%nblks_e,2) )
      IF (ALLOCATED(p_patch%edges%cell_blk)) DEALLOCATE(p_patch%edges%cell_blk)
      ALLOCATE( p_patch%edges%cell_blk(nproma,p_patch%nblks_e,2) )
      IF (ALLOCATED(p_patch%edges%vertex_idx)) DEALLOCATE(p_patch%edges%vertex_idx)
      ALLOCATE( p_patch%edges%vertex_idx(nproma,p_patch%nblks_e,4) )
      IF (ALLOCATED(p_patch%edges%vertex_blk)) DEALLOCATE(p_patch%edges%vertex_blk)
      ALLOCATE( p_patch%edges%vertex_blk(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%tangent_orientation(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%quad_idx(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%quad_blk(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%butterfly_idx(nproma,p_patch%nblks_e,2,2) )
      ALLOCATE( p_patch%edges%butterfly_blk(nproma,p_patch%nblks_e,2,2) )
      ALLOCATE( p_patch%edges%quad_orientation(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%center(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%primal_normal(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%dual_normal(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%primal_normal_cell(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%dual_normal_cell(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%primal_normal_vert(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%dual_normal_vert(nproma,p_patch%nblks_e,4) )
      ALLOCATE( p_patch%edges%primal_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%inv_primal_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%dual_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%inv_dual_edge_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%edge_vert_length(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%inv_vert_vert_length(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%edge_cell_length(nproma,p_patch%nblks_e,2) )
      ALLOCATE( p_patch%edges%area_edge(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%quad_area(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%f_e(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%fn_e(nproma,p_patch%nblks_e) )
      ALLOCATE( p_patch%edges%ft_e(nproma,p_patch%nblks_e) )
    ENDIF

    !
    ! !grid verts
    !
    IF (iopmode /= 2) THEN
      ALLOCATE( p_patch%verts%phys_id(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%neighbor_idx(nproma,p_patch%nblks_v,max_vertex_connectivity) )
      ALLOCATE( p_patch%verts%neighbor_blk(nproma,p_patch%nblks_v,max_vertex_connectivity) )
      IF (ALLOCATED(p_patch%verts%cell_idx)) DEALLOCATE(p_patch%verts%cell_idx)
      ALLOCATE( p_patch%verts%cell_idx(nproma,p_patch%nblks_v,max_vertex_connectivity) )
      IF (ALLOCATED(p_patch%verts%cell_blk)) DEALLOCATE(p_patch%verts%cell_blk)
      ALLOCATE( p_patch%verts%cell_blk(nproma,p_patch%nblks_v,max_vertex_connectivity) )
      ALLOCATE( p_patch%verts%edge_idx(nproma,p_patch%nblks_v,max_vertex_connectivity) )
      ALLOCATE( p_patch%verts%edge_blk(nproma,p_patch%nblks_v,max_vertex_connectivity) )
      ALLOCATE( p_patch%verts%edge_orientation(nproma,p_patch%nblks_v,max_vertex_connectivity) )
      IF (ALLOCATED(p_patch%verts%num_edges)) &
           DEALLOCATE(p_patch%verts%num_edges)
      ALLOCATE( p_patch%verts%num_edges(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%dual_area(nproma,p_patch%nblks_v) )
      ALLOCATE( p_patch%verts%f_v(nproma,p_patch%nblks_v) )
    ENDIF

    !
    ! ! decomp_info for cells, edges and verts
    !
    IF (iopmode /= 3) THEN

      CALL allocate_decomp_info(p_patch%cells%decomp_info, &
        &                  p_patch%n_patch_cells, &
        &                  p_patch%nblks_c)  ! this size is used for broadcasting,
                                             ! should not include the local ghost cell
       !  &                  p_patch%alloc_cell_blocks)
      CALL allocate_decomp_info(p_patch%edges%decomp_info, &
        &                  p_patch%n_patch_edges, &
        &                  p_patch%nblks_e)
      CALL allocate_decomp_info(p_patch%verts%decomp_info, &
        &                  p_patch%n_patch_verts, &
        &                  p_patch%nblks_v)

      CALL init_decomp_info(p_patch%cells%decomp_info, p_patch%npromz_c, &
        &                   p_patch%nblks_c, p_patch%n_patch_cells_g)
      CALL init_decomp_info(p_patch%edges%decomp_info, p_patch%npromz_e, &
        &                   p_patch%nblks_e, p_patch%n_patch_edges_g)
      CALL init_decomp_info(p_patch%verts%decomp_info, p_patch%npromz_v, &
        &                   p_patch%nblks_v, p_patch%n_patch_verts_g)
    END IF

    ! Set all newly allocated arrays to 0

    IF (iopmode /= 2) THEN
      p_patch%cells%phys_id = 0
      p_patch%cells%edge_orientation = 0._wp
      p_patch%cells%area = 0._wp
      p_patch%cells%f_c = 0._wp

      p_patch%edges%phys_id = 0
      p_patch%edges%cell_idx = 0
      p_patch%edges%cell_blk = 0
      p_patch%edges%vertex_idx = 0
      p_patch%edges%vertex_blk = 0
      p_patch%edges%tangent_orientation = 0._wp
      p_patch%edges%quad_idx = 0
      p_patch%edges%quad_blk = 0
      p_patch%edges%butterfly_idx = 0
      p_patch%edges%butterfly_blk = 0
      p_patch%edges%quad_orientation = 0._wp
      p_patch%edges%center(:,:)%lon = 0._wp
      p_patch%edges%center(:,:)%lat = 0._wp
      p_patch%edges%primal_normal(:,:)%v1 = 0._wp
      p_patch%edges%primal_normal(:,:)%v2 = 0._wp
      p_patch%edges%dual_normal(:,:)%v1 = 0._wp
      p_patch%edges%dual_normal(:,:)%v2 = 0._wp
      p_patch%edges%primal_normal_cell(:,:,:)%v1 = 0._wp
      p_patch%edges%primal_normal_cell(:,:,:)%v2 = 0._wp
      p_patch%edges%dual_normal_cell(:,:,:)%v1 = 0._wp
      p_patch%edges%dual_normal_cell(:,:,:)%v2 = 0._wp
      p_patch%edges%primal_normal_vert(:,:,:)%v1 = 0._wp
      p_patch%edges%primal_normal_vert(:,:,:)%v2 = 0._wp
      p_patch%edges%dual_normal_vert(:,:,:)%v1 = 0._wp
      p_patch%edges%dual_normal_vert(:,:,:)%v2 = 0._wp
      p_patch%edges%primal_edge_length = 0._wp
      p_patch%edges%inv_primal_edge_length = 0._wp
      p_patch%edges%dual_edge_length = 0._wp
      p_patch%edges%inv_dual_edge_length = 0._wp
      p_patch%edges%edge_vert_length = 0._wp
      p_patch%edges%inv_vert_vert_length = 0._wp
      p_patch%edges%edge_cell_length = 0._wp
      p_patch%edges%area_edge = 0._wp
      p_patch%edges%quad_area = 0._wp
      p_patch%edges%f_e = 0._wp
      p_patch%edges%fn_e = 0._wp
      p_patch%edges%ft_e = 0._wp

      p_patch%verts%phys_id = 0
      p_patch%verts%neighbor_idx = 0
      p_patch%verts%neighbor_blk = 0
      p_patch%verts%cell_idx = 0
      p_patch%verts%cell_blk = 0
      p_patch%verts%edge_idx = 0
      p_patch%verts%edge_blk = 0
      p_patch%verts%edge_orientation = 0._wp
      p_patch%verts%num_edges = 0
      p_patch%verts%dual_area = 0._wp
      p_patch%verts%f_v = 0._wp
    ENDIF

    IF (iopmode /= 2) CALL allocate_patch_cartesian( p_patch )

  CONTAINS

    SUBROUTINE allocate_decomp_info(decomp_info, n, n_blk)
      TYPE(t_grid_domain_decomp_info), INTENT(OUT) :: decomp_info
      INTEGER, INTENT(IN) :: n, n_blk

      ALLOCATE( decomp_info%decomp_domain(nproma,n_blk) )
      ! aliasing the halo_level to decomp_domain
      decomp_info%halo_level => decomp_info%decomp_domain
      ALLOCATE( decomp_info%owner_mask(nproma,n_blk) )
      ALLOCATE( decomp_info%glb_index(n) )
      ALLOCATE( decomp_info%owner_local(n))
    END SUBROUTINE allocate_decomp_info

    SUBROUTINE init_decomp_info(decomp_info, npromz, n_blk, n_g)
      TYPE(t_grid_domain_decomp_info), INTENT(INOUT) :: decomp_info
      INTEGER, INTENT(IN) :: npromz, n_blk, n_g

      INTEGER :: j

      ! Set values which are needed for parallel runs
      ! to the correct values for a single patch owner
      ! decomp_domain/owner mask:
      ! Everywhere 0 or .true. with the exception of unused entries

      decomp_info%decomp_domain = 0
      decomp_info%owner_mask = .TRUE.

      decomp_info%decomp_domain(npromz+1:nproma,n_blk) = -1
      decomp_info%owner_mask(npromz+1:nproma,n_blk) = .FALSE.

      ! The following arrays are currently never needed for non parallel runs,
      ! we set them nonetheless.
      ! MoHa: ...do they really need to be initialised?

      DO j = 1, SIZE(decomp_info%glb_index(:))
        decomp_info%glb_index(j) = j
      ENDDO
      decomp_info%owner_local(:) = 0

      CALL init_glb2loc_index_lookup(decomp_info%glb2loc_index, n_g)

    END SUBROUTINE init_decomp_info

  END SUBROUTINE allocate_remaining_patch
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !! Allocates complete patch by calling allocate_basic_patch/allocate_remaining_patch
  !!
  SUBROUTINE allocate_patch(p_patch)

    TYPE(t_patch), INTENT(inout) :: p_patch

    CALL allocate_basic_patch(p_patch)
    CALL allocate_remaining_patch(p_patch,1)

  END SUBROUTINE allocate_patch
  !-------------------------------------------------------------------------


END MODULE mo_alloc_patches


