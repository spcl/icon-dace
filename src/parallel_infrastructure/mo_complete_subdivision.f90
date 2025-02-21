! This module provides all routines for dividing patches
! (including interpolation state) and setting up communication.
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

MODULE mo_complete_subdivision
  !-------------------------------------------------------------------------
  USE mo_impl_constants,     ONLY: min_rlcell, min_rledge, &
    & min_rlcell_int, min_rledge_int, max_phys_dom
  USE mo_cdi_constants,      ONLY: grid_cell, grid_edge, grid_vertex
  USE mo_exception,          ONLY: finish, message, message_text

  USE mo_model_domain,       ONLY: t_patch, p_patch, p_patch_local_parent, &
    &                              p_phys_patch
  USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info, &
    &                              get_valid_local_index, &
    &                              t_glb2loc_index_lookup, &
    &                              partidx_of_elem_uniform_deco, &
    &                              uniform_partition
  USE mo_mpi,                ONLY: proc_split, p_sum
  USE mo_util_string,        ONLY: int2string






  USE mo_mpi,                ONLY: p_comm_work, my_process_is_mpi_test, &
    & my_process_is_mpi_seq, &
    & p_pe_work, p_n_work, p_bcast, &
    & p_comm_work_2_test, num_test_procs, p_comm_remote_size, &
    & p_isend, p_irecv, p_wait

  USE mo_communication,      ONLY: blk_no, idx_no, idx_1d, &
    &                              setup_comm_gather_pattern, t_comm_gather_pattern, &
    &                              ASSIGNMENT(=), delete_comm_gather_pattern, &
    &                              delete_comm_pattern, t_comm_pattern
  USE mo_communication_types, ONLY: xfer_list
  USE mo_communication_factory, ONLY: setup_comm_pattern, setup_comm_pattern2

  USE mo_impl_constants_grf, ONLY: grf_bdyintp_start_c, grf_bdyintp_start_e,  &
    & grf_fbk_start_c, grf_fbk_start_e, grf_bdywidth_c, &
    & grf_bdywidth_e
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, n_phys_dom
  USE mo_dist_dir,            ONLY: dist_dir_get_owners
  USE ppm_distributed_array, ONLY: global_array_desc, dist_mult_array, &



       sync_mode_passive_target, &

       dist_mult_array_new, dist_mult_array_delete, &
       dist_mult_array_get, dist_mult_array_local_ptr, &
       dist_mult_array_expose
  USE mo_kind, ONLY: wp
  USE mo_sync, ONLY: sync_c, sync_patch_array
  IMPLICIT NONE

  PRIVATE

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: copy_processor_splitting
  PUBLIC :: finalize_decomposition
  PUBLIC :: setup_phys_patches
  PUBLIC :: complete_parallel_setup
  PUBLIC :: generate_comm_pat_cvec1




CONTAINS

  !-------------------------------------------------------------------------
  !>
  !!  In case of a test run: Copies processor splitting to test PE
  SUBROUTINE copy_processor_splitting(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)


  END SUBROUTINE copy_processor_splitting

  !-------------------------------------------------------------------------
  !>
  !!  Sets the communicators in the patches if these have been read.
  SUBROUTINE set_patch_communicators(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch(n_dom_start:)

    INTEGER jc, jgc, jg, jgp, n_proc_total, comm, mpierr
    INTEGER, ALLOCATABLE :: patch_no(:)

    IF (my_process_is_mpi_seq()) THEN
      p_patch(:)%comm = p_comm_work
      RETURN
    ENDIF


  END SUBROUTINE set_patch_communicators

  !-----------------------------------------------------------------------------
  ! sets communication patterns and parent-child relationships.
  !
  SUBROUTINE generate_comm_pat_cvec1(patch, is_ocean_decomposition)

    TYPE(t_patch), INTENT(INOUT) :: patch(n_dom_start:)
    LOGICAL, INTENT(IN) :: is_ocean_decomposition

    INTEGER :: jg

    IF (is_ocean_decomposition .AND. (n_dom > n_dom_start)) &
      CALL finish('generate_comm_pat_cvec1', &
        &         'functionality with local parent patch not implemented')

    DO jg = n_dom_start, n_dom

      ! Set communication patterns for boundary exchange
      CALL set_comm_pat_bound_exch(patch(jg))

      IF (jg > n_dom_start) THEN

        CALL set_comm_pat_bound_exch(p_patch_local_parent(jg))
      ENDIF

    ENDDO

  END SUBROUTINE generate_comm_pat_cvec1

  !-----------------------------------------------------------------------------
  ! sets communication patterns and parent-child relationships.
  !
  SUBROUTINE complete_parallel_setup(patch, is_ocean_decomposition)

    TYPE(t_patch), INTENT(INOUT) :: patch(n_dom_start:)
    LOGICAL, INTENT(IN) :: is_ocean_decomposition

    INTEGER :: jg, jgp, i

    IF (is_ocean_decomposition .AND. (n_dom > n_dom_start)) &
      CALL finish('complete_parallel_setup', &
        &         'functionality with local parent patch not implemented')

    DO jg = n_dom_start, n_dom

      jgp = patch(jg)%parent_id

      ! Rebuild communication patterns for boundary exchange
      CALL delete_comm_pattern(patch(jg)%comm_pat_c)
      CALL delete_comm_pattern(patch(jg)%comm_pat_v)
      CALL delete_comm_pattern(patch(jg)%comm_pat_e)
      CALL delete_comm_pattern(patch(jg)%comm_pat_c1)
      CALL set_comm_pat_bound_exch(patch(jg))

      ! Set communication patterns for gathering on proc 0
      CALL set_comm_pat_gather(patch(jg))

      ! Create communication pattern with test processes
      IF (jg > n_dom_start) THEN

        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_c)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_v)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_e)
        CALL delete_comm_pattern(p_patch_local_parent(jg)%comm_pat_c1)
        CALL set_comm_pat_bound_exch(p_patch_local_parent(jg))

        CALL set_glb_loc_comm(patch(jgp), p_patch_local_parent(jg), &
          &                   patch(jg)%parent_child_index)
      ENDIF

    ENDDO

    CALL set_patch_communicators(patch)

  END SUBROUTINE complete_parallel_setup

  !-----------------------------------------------------------------------------

  SUBROUTINE finalize_decomposition(patch, is_ocean_decomposition)

    TYPE(t_patch), INTENT(INOUT) :: patch(n_dom_start:)
    LOGICAL, INTENT(IN) :: is_ocean_decomposition

    INTEGER :: jg

    IF (is_ocean_decomposition .AND. (n_dom > n_dom_start)) &
      CALL finish('finalize_decomposition_oce', &
        &         'functionality with local parent patch not implemented')

    ! Remap indices in patches and local parents

    DO jg = n_dom_start, n_dom

      CALL remap_patch_indices(patch(jg))

      IF(jg>n_dom_start) THEN
        CALL remap_patch_indices(p_patch_local_parent(jg))
      ENDIF

    ENDDO

  END SUBROUTINE finalize_decomposition

  !-------------------------------------------------------------------------------------------------
  !>
  !! Remaps negative index entries to valid ones
  !!

  SUBROUTINE remap_patch_indices(p_patch)

    TYPE(t_patch), INTENT(INOUT) :: p_patch

    INTEGER :: i, j, jb, jl, ilc1, ibc1, ilc2, ibc2, ilc3, ibc3, irl0, irl1, irl2, irl3

    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,p_patch%cells%max_connectivity

!CDIR IEXPAND
        CALL remap_index(p_patch%cells%decomp_info%glb2loc_index, &
          & p_patch%cells%neighbor_idx(jl,jb,i),    &
          & p_patch%cells%neighbor_blk(jl,jb,i))

        ! edge_idx and vertex_idx should not need a remap !!!
        ! This is only left here if there should change something in the decomposition
!CDIR IEXPAND
        CALL remap_index(p_patch%edges%decomp_info%glb2loc_index, &
          & p_patch%cells%edge_idx(jl,jb,i),        &
          & p_patch%cells%edge_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%verts%decomp_info%glb2loc_index, &
          & p_patch%cells%vertex_idx(jl,jb,i),      &
          & p_patch%cells%vertex_blk(jl,jb,i))

      ENDDO
    ENDDO

    ! ensure that cells%neighbor_idx lies in the correct grid row along the lateral boundary
    DO j = 1, p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (p_patch%cells%refin_ctrl(jl,jb) > 0 .AND.             &
          p_patch%cells%refin_ctrl(jl,jb) <= grf_bdywidth_c) THEN

        ilc1 = p_patch%cells%neighbor_idx(jl,jb,1)
        ibc1 = p_patch%cells%neighbor_blk(jl,jb,1)
        ilc2 = p_patch%cells%neighbor_idx(jl,jb,2)
        ibc2 = p_patch%cells%neighbor_blk(jl,jb,2)
        ilc3 = p_patch%cells%neighbor_idx(jl,jb,3)
        ibc3 = p_patch%cells%neighbor_blk(jl,jb,3)

        irl0 = p_patch%cells%refin_ctrl(jl,jb)

        IF (ilc1 > 0 .AND. ibc1 > 0) THEN
          irl1 = p_patch%cells%refin_ctrl(ilc1,ibc1)
          IF ( irl1 > 0 .AND. ABS(irl0 - irl1) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,1) = jl
            p_patch%cells%neighbor_blk(jl,jb,1) = jb
          ENDIF
        ENDIF

        IF (ilc2 > 0 .AND. ibc2 > 0) THEN
          irl2 = p_patch%cells%refin_ctrl(ilc2,ibc2)
          IF ( irl2 > 0 .AND. ABS(irl0 - irl2) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,2) = jl
            p_patch%cells%neighbor_blk(jl,jb,2) = jb
          ENDIF
        ENDIF

        IF (ilc3 > 0 .AND. ibc3 > 0) THEN
          irl3 = p_patch%cells%refin_ctrl(ilc3,ibc3)
          IF ( irl3 > 0 .AND. ABS(irl0 - irl3) >1 ) THEN
            p_patch%cells%neighbor_idx(jl,jb,3) = jl
            p_patch%cells%neighbor_blk(jl,jb,3) = jb
          ENDIF
        ENDIF

      ENDIF
    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (p_patch%edges%refin_ctrl(jl,jb) /= 1) THEN
        DO i=1,2
!CDIR IEXPAND
          CALL remap_index(p_patch%cells%decomp_info%glb2loc_index, &
            & p_patch%edges%cell_idx(jl,jb,i),        &
            & p_patch%edges%cell_blk(jl,jb,i))
        ENDDO
      ENDIF

      DO i=1,4
!CDIR IEXPAND
        CALL remap_index(p_patch%verts%decomp_info%glb2loc_index, &
          & p_patch%edges%vertex_idx(jl,jb,i),      &
          & p_patch%edges%vertex_blk(jl,jb,i))
      ENDDO

      DO i=1,4
!CDIR IEXPAND
        CALL remap_index(p_patch%edges%decomp_info%glb2loc_index, &
          & p_patch%edges%quad_idx(jl,jb,i),        &
          & p_patch%edges%quad_blk(jl,jb,i))
      ENDDO

      ! ensure that edges%cell_idx lies in the correct grid row along the lateral boundary
      IF (p_patch%edges%refin_ctrl(jl,jb) >= 2 .AND.            &
          p_patch%edges%refin_ctrl(jl,jb) <= grf_bdywidth_e) THEN

        ilc1 = p_patch%edges%cell_idx(jl,jb,1)
        ibc1 = p_patch%edges%cell_blk(jl,jb,1)
        ilc2 = p_patch%edges%cell_idx(jl,jb,2)
        ibc2 = p_patch%edges%cell_blk(jl,jb,2)
        irl1 = p_patch%cells%refin_ctrl(ilc1,ibc1)
        irl2 = p_patch%cells%refin_ctrl(ilc2,ibc2)

        IF (MOD(p_patch%edges%refin_ctrl(jl,jb),2)==0) THEN
          IF (irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            p_patch%edges%cell_idx(jl,jb,1) = p_patch%edges%cell_idx(jl,jb,2)
            p_patch%edges%cell_blk(jl,jb,1) = p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2) THEN
            p_patch%edges%cell_idx(jl,jb,2) = p_patch%edges%cell_idx(jl,jb,1)
            p_patch%edges%cell_blk(jl,jb,2) = p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ELSE
          IF (irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl1 /= p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            p_patch%edges%cell_idx(jl,jb,1) = p_patch%edges%cell_idx(jl,jb,2)
            p_patch%edges%cell_blk(jl,jb,1) = p_patch%edges%cell_blk(jl,jb,2)
          ENDIF
          IF (irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2 .AND. &
              irl2 /= p_patch%edges%refin_ctrl(jl,jb)/2+1) THEN
            p_patch%edges%cell_idx(jl,jb,2) = p_patch%edges%cell_idx(jl,jb,1)
            p_patch%edges%cell_blk(jl,jb,2) = p_patch%edges%cell_blk(jl,jb,1)
          ENDIF
        ENDIF

      ENDIF

    ENDDO

    !---------------------------------------------------------------------------------------

    DO j = 1,p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,p_patch%verts%max_connectivity

!CDIR IEXPAND
        CALL remap_index(p_patch%verts%decomp_info%glb2loc_index, &
          & p_patch%verts%neighbor_idx(jl,jb,i),    &
          & p_patch%verts%neighbor_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%edges%decomp_info%glb2loc_index, &
          & p_patch%verts%edge_idx(jl,jb,i),        &
          & p_patch%verts%edge_blk(jl,jb,i))

!CDIR IEXPAND
        CALL remap_index(p_patch%cells%decomp_info%glb2loc_index, &
          & p_patch%verts%cell_idx(jl,jb,i),        &
          & p_patch%verts%cell_blk(jl,jb,i))
      ENDDO

    ENDDO

  END SUBROUTINE remap_patch_indices

  !-------------------------------------------------------------------------------------------------
  !> Sets the communication patterns for boundary exchange

  SUBROUTINE set_comm_pat_bound_exch(p)

    TYPE(t_patch), INTENT(INOUT):: p

    INTEGER, ALLOCATABLE :: owner_c(:), owner_e(:), owner_v(:)
    INTEGER :: jc

    ALLOCATE(owner_c(p%n_patch_cells))
    ALLOCATE(owner_e(p%n_patch_edges))
    ALLOCATE(owner_v(p%n_patch_verts))

    ! Set the owner arrays for cells/edges/verts which have to be transferred.
    ! The following setting always transfers edges/verts if they are owned
    ! (according to the MIN/MAX PE setting above) by another PE, which implies
    ! that boundary edges/verts (with flag2_e/flag2_v = 1) can be excluded from
    ! prognostic computations if they are followed by a synchronization call.

    owner_c(:) = MERGE(p%cells%decomp_info%owner_local(:), -1, &
      &                p%cells%decomp_info%owner_local(:) /= p_pe_work)

    owner_e(:) = MERGE(p%edges%decomp_info%owner_local(:), -1, &
      &                p%edges%decomp_info%owner_local(:) /= p_pe_work)

    owner_v(:) = MERGE(p%verts%decomp_info%owner_local(:), -1, &
      &                p%verts%decomp_info%owner_local(:) /= p_pe_work)

    ! Set communication patterns for boundary exchange
    CALL setup_comm_pattern(p%n_patch_cells, owner_c, p%cells%decomp_info%glb_index, &
      & p%cells%decomp_info%glb2loc_index, p%n_patch_cells, &
      & p%cells%decomp_info%owner_local, p%cells%decomp_info%glb_index, &
      & p%comm_pat_c, inplace=.TRUE.)

    CALL setup_comm_pattern(p%n_patch_edges, owner_e, p%edges%decomp_info%glb_index, &
      & p%edges%decomp_info%glb2loc_index, p%n_patch_edges, &
      & p%edges%decomp_info%owner_local, p%edges%decomp_info%glb_index, &
      p%comm_pat_e, inplace=.TRUE.)

    CALL setup_comm_pattern(p%n_patch_verts, owner_v, p%verts%decomp_info%glb_index, &
      & p%verts%decomp_info%glb2loc_index, p%n_patch_verts, &
      & p%verts%decomp_info%owner_local, p%verts%decomp_info%glb_index, &
      & p%comm_pat_v, inplace=.TRUE.)

    DEALLOCATE(owner_e, owner_v)

    ! Set reduced communication pattern containing only level-1 halo cells (immediate neighbors)
    jc = idx_1d(p%cells%end_idx(min_rlcell_int-1,1), &
                p%cells%end_blk(min_rlcell_int-1,1))
    CALL setup_comm_pattern(jc, owner_c(1:jc), p%cells%decomp_info%glb_index, &
      & p%cells%decomp_info%glb2loc_index, p%n_patch_cells, &
      & p%cells%decomp_info%owner_local, p%cells%decomp_info%glb_index, &
      & p%comm_pat_c1, inplace=.TRUE.)

    DEALLOCATE(owner_c)

  END SUBROUTINE set_comm_pat_bound_exch

  !-------------------------------------------------------------------------------------------------
  !> Sets the gather communication patterns of a patch

  SUBROUTINE set_comm_pat_gather(p)

    TYPE(t_patch), INTENT(INOUT):: p

    CALL setup_comm_gather_pattern( &
      p%cells%decomp_info%glb2loc_index%global_size, &
      p%cells%decomp_info%owner_local, &
      p%cells%decomp_info%glb_index, p%comm_pat_gather_c)
    CALL setup_comm_gather_pattern( &
      p%verts%decomp_info%glb2loc_index%global_size, &
      p%verts%decomp_info%owner_local, &
      p%verts%decomp_info%glb_index, p%comm_pat_gather_v)
    CALL setup_comm_gather_pattern( &
      p%edges%decomp_info%glb2loc_index%global_size, &
      p%edges%decomp_info%owner_local, &
      p%edges%decomp_info%glb_index, p%comm_pat_gather_e)

  END SUBROUTINE set_comm_pat_gather

  !-------------------------------------------------------------------------------------------------
  !
  !> Sets up communication patterns between global and local parent patches

  SUBROUTINE set_glb_loc_comm(p_pglb, p_ploc, i_chidx)
    TYPE(t_patch), INTENT(IN)    :: p_pglb !> global parent
    TYPE(t_patch), INTENT(INOUT) :: p_ploc !> local parent
    INTEGER, INTENT(IN) :: i_chidx

    INTEGER, ALLOCATABLE :: owner(:)
    LOGICAL, ALLOCATABLE :: mask(:)
    INTEGER :: j, je, icid, jb, jl, max_size

    ! Please note:
    ! For creating communication patterns for different amount of data to be transferred
    ! (e.g. only the boundary interpolation zone), create a new communication pattern
    ! by copying the code below and adjusting the limits in the calculation of js/je.

    !-----------------------------------------------------------------------------------------------

    ! child ID
    icid = p_pglb%child_id(i_chidx)

    ! Communication global -> local
    ! Only one pattern is set which can be used everywhere since it doesn't matter
    ! if more cells/edges than needed are set in the local parent.
    ! Please note that only the inner area of the local patch is set and not the boundary cells
    ! (which are needed only for gradient calculation in the moment).

    ! ... cells

    max_size = MAX(p_ploc%n_patch_cells, p_ploc%n_patch_edges, &
      &            p_pglb%n_patch_cells, p_pglb%n_patch_edges)
    ALLOCATE(owner(max_size), mask(max_size))

    je = idx_1d(p_ploc%cells%end_idx(min_rlcell_int,i_chidx), &
      &         p_ploc%cells%end_blk(min_rlcell_int,i_chidx))

    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_ploc%cells%child_id(jl,jb)   == icid .AND. &
        &       p_ploc%cells%refin_ctrl(jl,jb) <= grf_bdyintp_start_c
    ENDDO
    mask(je+1:p_ploc%n_patch_cells) = .FALSE.

    owner(1:p_ploc%n_patch_cells) = &
         dist_dir_get_owners(p_pglb%cells%decomp_info%owner_dist_dir, &
         p_ploc%cells%decomp_info%glb_index(:), &
         mask(1:p_ploc%n_patch_cells))

    CALL setup_comm_pattern(p_ploc%n_patch_cells, owner(1:p_ploc%n_patch_cells), &
      p_ploc%cells%decomp_info%glb_index, &
      p_pglb%cells%decomp_info%glb2loc_index, &
      p_pglb%n_patch_cells, &
      p_pglb%cells%decomp_info%owner_local, &
      p_pglb%cells%decomp_info%glb_index, &
      p_ploc%comm_pat_glb_to_loc_c)

    ! ... edges

    je = idx_1d(p_ploc%edges%end_idx(min_rledge_int,i_chidx), &
      &         p_ploc%edges%end_blk(min_rledge_int,i_chidx))

    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_ploc%edges%child_id(jl,jb) == icid .AND. &
        &       p_ploc%edges%refin_ctrl(jl,jb) <= grf_bdyintp_start_e
    ENDDO
    mask(je+1:p_ploc%n_patch_edges) = .FALSE.

    owner(1:p_ploc%n_patch_edges) = &
         dist_dir_get_owners(p_pglb%edges%decomp_info%owner_dist_dir, &
         p_ploc%edges%decomp_info%glb_index(:), &
         mask(1:p_ploc%n_patch_edges))

    CALL setup_comm_pattern(p_ploc%n_patch_edges, owner(1:p_ploc%n_patch_edges), &
      p_ploc%edges%decomp_info%glb_index,  &
      p_pglb%edges%decomp_info%glb2loc_index, &
      p_pglb%n_patch_edges, &
      p_pglb%edges%decomp_info%owner_local, &
      p_pglb%edges%decomp_info%glb_index, &
      p_ploc%comm_pat_glb_to_loc_e)

    !-----------------------------------------------------------------------------------------------

    ! Communication local -> global
    ! Here it might get necessary to have different patterns for different start levels
    ! of the copy to the global parent since we may not overwrite arbitrary values there.
    ! Currently only one pattern for feedback is needed (starting at grf_fbk_start_c/e).

    ! ... cells

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      je = idx_1d(p_pglb%cells%end_idx(min_rlcell,p_pglb%n_childdom), &
        &         p_pglb%cells%end_blk(min_rlcell,p_pglb%n_childdom))
    ELSE
      je = idx_1d(p_pglb%cells%end_idx(min_rlcell_int,p_pglb%n_childdom), &
        &         p_pglb%cells%end_blk(min_rlcell_int,p_pglb%n_childdom))
    ENDIF

    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_pglb%cells%child_id(jl,jb) == icid .AND. &
        &       p_pglb%cells%refin_ctrl(jl,jb) <= grf_fbk_start_c .AND. &
        &       p_pglb%cells%refin_ctrl(jl,jb) >= min_rlcell_int
    ENDDO
    mask(je+1:p_pglb%n_patch_cells) = .FALSE.

    owner(1:p_pglb%n_patch_cells) = &
         dist_dir_get_owners(p_ploc%cells%decomp_info%owner_dist_dir, &
         p_pglb%cells%decomp_info%glb_index(:), &
         mask(1:p_pglb%n_patch_cells))

    CALL setup_comm_pattern(p_pglb%n_patch_cells, owner(1:p_pglb%n_patch_cells), &
      p_pglb%cells%decomp_info%glb_index, &
      p_ploc%cells%decomp_info%glb2loc_index, &
      p_ploc%n_patch_cells, &
      p_ploc%cells%decomp_info%owner_local, &
      p_ploc%cells%decomp_info%glb_index, &
      p_ploc%comm_pat_loc_to_glb_c_fbk)

    ! ... edges

    IF (p_pglb%id > 0) THEN  ! include halo points belonging to nest overlap points
      je = idx_1d(p_pglb%edges%end_idx(min_rledge,i_chidx), &
        &         p_pglb%edges%end_blk(min_rledge,i_chidx))
    ELSE
      je = idx_1d(p_pglb%edges%end_idx(min_rledge_int,i_chidx), &
        &         p_pglb%edges%end_blk(min_rledge_int,i_chidx))
    ENDIF

    DO j = 1, je
      jb = blk_no(j) ! Block index
      jl = idx_no(j) ! Line  index
      mask(j) = p_pglb%edges%child_id(jl,jb) == icid .AND. &
        &       p_pglb%edges%refin_ctrl(jl,jb) <= grf_fbk_start_e .AND. &
        &       p_pglb%edges%refin_ctrl(jl,jb) >= min_rledge_int
    ENDDO
    mask(je+1:p_pglb%n_patch_edges) = .FALSE.

    owner(1:p_pglb%n_patch_edges) = &
         dist_dir_get_owners(p_ploc%edges%decomp_info%owner_dist_dir, &
         p_pglb%edges%decomp_info%glb_index(:), &
         mask(1:p_pglb%n_patch_edges))

    CALL setup_comm_pattern(p_pglb%n_patch_edges, &
      owner(1:p_pglb%n_patch_edges), &
      p_pglb%edges%decomp_info%glb_index, &
      p_ploc%edges%decomp_info%glb2loc_index, &
      p_ploc%n_patch_edges, &
      p_ploc%edges%decomp_info%owner_local, &
      p_ploc%edges%decomp_info%glb_index, &
      p_ploc%comm_pat_loc_to_glb_e_fbk)

    DEALLOCATE(owner, mask)
  END SUBROUTINE set_glb_loc_comm

  !-------------------------------------------------------------------------
  !>
  !! This routine sets up a the physical patches
  !! It works on the divided patch state, so it can be used after a restore also

  SUBROUTINE setup_phys_patches
    INTEGER :: jp, jg, i
    CHARACTER(LEN=*), PARAMETER :: routine = 'setup_phys_patches'

    INTEGER :: n_patch_cells(max_phys_dom), n_patch_verts(max_phys_dom), &
      &        n_patch_edges(max_phys_dom)
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_cells(max_phys_dom), &
      &                            comm_pat_gather_verts(max_phys_dom), &
      &                            comm_pat_gather_edges(max_phys_dom)

    p_phys_patch(:)%logical_id = -1
    n_patch_cells(:) = 0
    n_patch_verts(:) = 0
    n_patch_edges(:) = 0

    DO jg = 1, n_dom
      CALL setup_phys_patches_cve("cells, jg="//int2string(jg), &
        &                         p_patch(jg)%n_patch_cells_g, &
        &                         p_patch(jg)%n_patch_cells, &
        &                         p_patch(jg)%cells%decomp_info, &
        &                         p_patch(jg)%cells%phys_id, jg, &
        &                         .TRUE., n_patch_cells(:), &
        &                         comm_pat_gather_cells(:))
      CALL setup_phys_patches_cve("edges, jg="//int2string(jg), &
        &                         p_patch(jg)%n_patch_verts_g, &
        &                         p_patch(jg)%n_patch_verts, &
        &                         p_patch(jg)%verts%decomp_info, &
        &                         p_patch(jg)%verts%phys_id, jg, &
        &                         .FALSE., n_patch_verts(:), &
        &                         comm_pat_gather_verts(:))
      CALL setup_phys_patches_cve("verts, jg="//int2string(jg), &
        &                         p_patch(jg)%n_patch_edges_g, &
        &                         p_patch(jg)%n_patch_edges, &
        &                         p_patch(jg)%edges%decomp_info, &
        &                         p_patch(jg)%edges%phys_id, jg, &
        &                         .FALSE., n_patch_edges(:), &
        &                         comm_pat_gather_edges(:))
    ENDDO

    ! MoHa: remark:
    ! We cannot pass p_phys_patch(:)%n_patch_cells and
    ! p_phys_patch(:)%comm_pat_gather_c directly to setup_phys_patches_cve
    ! because some compiler(intel and pgi) have problems with this type
    ! of array argument... Therefore, the temporary arrays n_patch_*(:) and
    ! comm_pat_gather_*(:) had to be used.
    DO i = 1, max_phys_dom
      p_phys_patch(i)%n_patch_cells = n_patch_cells(i)
      p_phys_patch(i)%comm_pat_gather_c = comm_pat_gather_cells(i)
      CALL delete_comm_gather_pattern(comm_pat_gather_cells(i))
      p_phys_patch(i)%n_patch_verts = n_patch_verts(i)
      p_phys_patch(i)%comm_pat_gather_v = comm_pat_gather_verts(i)
      CALL delete_comm_gather_pattern(comm_pat_gather_verts(i))
      p_phys_patch(i)%n_patch_edges = n_patch_edges(i)
      p_phys_patch(i)%comm_pat_gather_e = comm_pat_gather_edges(i)
      CALL delete_comm_gather_pattern(comm_pat_gather_edges(i))
    ENDDO

    ! Get total number of physical patches

    DO jp = max_phys_dom, 1, -1
      IF(p_phys_patch(jp)%logical_id > 0) EXIT
    ENDDO

    n_phys_dom = jp

    ! Print info and check if there are no unused physical patches

    DO jp = 1, n_phys_dom
      WRITE(message_text,'(a,i4,a,i4)') &
        & 'Physical domain ',jp,' belongs to logical domain ',p_phys_patch(jp)%logical_id
      CALL message (routine, message_text)
      IF(p_phys_patch(jp)%logical_id < 0) THEN
        WRITE(message_text,'(a,i4)') 'Missing physical domain # ',jp
        CALL finish (routine, message_text)
      ENDIF
    ENDDO

  END SUBROUTINE setup_phys_patches

  SUBROUTINE setup_phys_patches_cve(description_str, n_g, n, decomp_info, phys_id, &
    &                               curr_patch_idx, set_logical_id, &
    &                               n_patch_cve, comm_pat_gather)
    CHARACTER(LEN=*), INTENT(IN) :: description_str
    INTEGER, INTENT(IN) :: n_g, n
    TYPE(t_grid_domain_decomp_info), INTENT(IN) :: decomp_info
    INTEGER, INTENT(IN) :: phys_id(:,:)
    INTEGER, INTENT(IN) :: curr_patch_idx
    LOGICAL, INTENT(IN) :: set_logical_id
    INTEGER, INTENT(INOUT) :: n_patch_cve(:)
    TYPE(t_comm_gather_pattern), INTENT(INOUT) :: comm_pat_gather(:)

    INTEGER, ALLOCATABLE :: owner_local(:)
    INTEGER :: temp_n_patch_cve(max_phys_dom)
    INTEGER :: i, ip

    ALLOCATE(owner_local(n))

    temp_n_patch_cve(:) = 0

    ! Fill with own values of phys_id

    DO i = 1, n
      ip = phys_id(idx_no(i),blk_no(i))
      IF (ip < 1 .OR. ip > max_phys_dom) &
        CALL finish("setup_phys_patches_cve, "//TRIM(description_str), "invalid phys_id "//int2string(ip))
      IF (decomp_info%owner_local(i) == p_pe_work) &
        temp_n_patch_cve(ip) = temp_n_patch_cve(ip) + 1
    END DO

    temp_n_patch_cve(:) = p_sum(temp_n_patch_cve(:), p_comm_work)

    IF (set_logical_id) THEN

      ! Check if no other patch uses the same phys_id
      IF (ANY((p_phys_patch(1:max_phys_dom)%logical_id /= -1) .AND. &
        &     (temp_n_patch_cve(:) /= 0))) &
        CALL finish("setup_phys_patches_cve, "//TRIM(description_str), &
          &         "invalid value in p_phys_patch(:)%logical_id")

      WHERE(temp_n_patch_cve(:) /= 0) &
        p_phys_patch(1:max_phys_dom)%logical_id = curr_patch_idx
    ELSE
      ! Check if no other patch uses the same phys_id
      IF (ANY((p_phys_patch(1:max_phys_dom)%logical_id /= curr_patch_idx) &
        &     .AND. (temp_n_patch_cve(:) /= 0))) &
        CALL finish("setup_phys_patches_cve, "//TRIM(description_str), &
          &         "invalid value in p_phys_patch(:)%logical_id")
    END IF

    DO ip = 1, max_phys_dom

      IF (p_phys_patch(ip)%logical_id /= curr_patch_idx) CYCLE

      n_patch_cve(ip) = temp_n_patch_cve(ip)

      DO i = 1, n
        owner_local(i) = &
          MERGE(decomp_info%owner_local(i), -1, &
            &   phys_id(idx_no(i),blk_no(i)) == ip)
      END DO

      CALL setup_comm_gather_pattern(n_g, owner_local(:), &
        &                            decomp_info%glb_index(:), &
        &                            comm_pat_gather(ip), .TRUE.)
    END DO

    DEALLOCATE(owner_local)

  END SUBROUTINE setup_phys_patches_cve

  !-------------------------------------------------------------------------
  !! Calculates local line/block indices l_idx, l_blk
  !! from global line/block indices g_idx, g_blk
  !! using the mapping in decomp_info
  !!

  !-------------------------------------------------------------------------
  !! Maps indices which point outside domain (returned as <0 by get_local_idx_blk)
  !! to valid ones which are from the same set of neighbors for the given
  !! cell/edge/vert
  !!
  SUBROUTINE remap_index(glb2loc_index, l_idx, l_blk)

    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: glb2loc_index
    INTEGER, INTENT(INOUT) :: l_idx, l_blk

    INTEGER :: j_l

    IF(l_idx>=0) RETURN ! Nothing to do

    j_l = MAX(get_valid_local_index(glb2loc_index, &
      &                             idx_1d(-l_idx, l_blk)), 1)
    l_idx = idx_no(j_l)
    l_blk = blk_no(j_l)

  END SUBROUTINE remap_index

  !-----------------------------------------------------------------------------
END MODULE mo_complete_subdivision

