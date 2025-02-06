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

MODULE mo_build_decomposition

  USE mo_kind,                     ONLY: wp
  USE mo_complete_subdivision, ONLY: finalize_decomposition,                  &
    &                                copy_processor_splitting,                &
    &                                complete_parallel_setup
  USE mo_setup_subdivision,    ONLY: decompose_domain
  USE mo_sync,                 ONLY: disable_sync_checks, enable_sync_checks, &
    &                                decomposition_statistics
  USE mo_grid_config,          ONLY: n_dom, n_dom_start
  USE mo_mpi,                  ONLY: my_process_is_mpi_parallel, p_pe_work,   &
    &                                p_comm_work, p_comm_work_test
  USE mo_loopindices,          ONLY: get_indices_e
  USE mo_model_domain,         ONLY: p_patch, t_pre_patch, t_patch_3d,        &
    &                                t_patch, p_patch_local_parent
  USE mo_reshuffle,            ONLY: reshuffle
  USE mo_run_config,           ONLY: msg_level
  USE mo_model_domimp_patches, ONLY: reorder_patch_refin_ctrl,                &
    &                                import_pre_patches, complete_patches, set_parent_loc_idx
  USE mo_parallel_config,      ONLY: p_test_run, l_test_openmp, num_io_procs, division_method
  USE mo_impl_constants,       ONLY: success, max_dom, min_rlcell_int
  USE mo_exception,            ONLY: finish, message, message_text
  USE mo_communication,        ONLY: blk_no, idx_no, idx_1d
  USE mo_util_string,          ONLY: int2string

  IMPLICIT NONE
  
  PUBLIC :: build_decomposition

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_build_decomposition'
  
CONTAINS

  !> Main routine for creating the domain decomposition (together with
  !> communication patterns etc.)
  !
  !  @note This routine is called for both: The ocean model and the
  !        atmo_model.
  !
  SUBROUTINE build_decomposition(num_lev,nshift,&
    &                            is_ocean_decomposition, patch_3d)
    
    INTEGER, INTENT(in)                 :: num_lev(max_dom),                &
      &                                    nshift(max_dom)
    LOGICAL, INTENT(in)                 :: is_ocean_decomposition
    TYPE(t_patch_3d), POINTER, OPTIONAL :: patch_3d
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//':build_decomposition'
    TYPE(t_pre_patch), ALLOCATABLE :: p_patch_pre(:)
    INTEGER                        :: error_status, jg, jgp
    !> If .true., read fields related to grid refinement from separate  grid files:
    LOGICAL                        :: lsep_grfinfo
    
    ! Check p_patch allocation status
    
    IF ( ALLOCATED(p_patch)) THEN
      CALL finish(TRIM(routine), 'p_patch already allocated')
    END IF
    
    ! Allocate p_patch array to start p_patch construction.
    !
    ! At the same time, we allocate the "p_patch_local_parent" which
    ! is the local portion of each p_patch's parent grid.
    ALLOCATE(p_patch             (n_dom_start  :n_dom), &
      & p_patch_local_parent(n_dom_start+1:n_dom), &
      & stat=error_status)
    IF (error_status/=success) CALL finish(routine, 'allocation of p_patch failed')
    
    ! --------------------------
    ! Work PEs subdivide patches
    ! --------------------------
                  
    ! compute domain decomposition on-the-fly        
    ALLOCATE(p_patch_pre(n_dom_start:n_dom))
    CALL import_pre_patches(p_patch_pre,num_lev,nshift,lsep_grfinfo)
    ! use internal domain decomposition algorithm
    CALL decompose_domain(p_patch, p_patch_pre)
    DEALLOCATE(p_patch_pre)

    ! compute the child/edge indices on parent patch from the
    ! corresponding information on parent indices stored in the child
    ! patch
    DO jg = n_dom_start,n_dom
      ! set child indices of patch:
      jgp = p_patch(jg)%parent_id
      IF (jgp >= n_dom_start) THEN
        CALL set_child_indices("patch "//TRIM(int2string(jgp))//" -> "//int2string(jg), &
          &                    p_patch(jg), p_patch(jgp), is_local_parent=.FALSE.)
        
        ! set child indices of local parent patch, but don't modify the
        ! "pc_idx" fields of the child patch:
        CALL set_child_indices("loc par patch "//TRIM(int2string(jgp))//" -> "//int2string(jg), &
          &                    p_patch(jg), p_patch_local_parent(jg), is_local_parent=.TRUE.)
      END IF
    END DO

    ! Complete information which is not yet read or calculated
    CALL complete_patches( p_patch, is_ocean_decomposition, lsep_grfinfo)

    ! In case of a test run: Copy processor splitting to test PE
    IF(p_test_run) CALL copy_processor_splitting(p_patch)
    !--------------------------------------------------------------------------------

    CALL finalize_decomposition(p_patch, is_ocean_decomposition)

    ! reorder patch_local_parents according to their refin_ctrl flags
    DO jg = n_dom_start+1, n_dom
      CALL reorder_patch_refin_ctrl(p_patch_local_parent(jg), p_patch(jg))
    ENDDO

    ! computes communication patterns (done after reordering)
    CALL complete_parallel_setup(p_patch, is_ocean_decomposition)

    IF(.NOT.p_test_run .AND. my_process_is_mpi_parallel()) THEN ! the call below hangs in test mode
      ! Print diagnostic information about domain decomposition
      DO jg = 1, n_dom
        CALL decomposition_statistics(p_patch(jg))
      ENDDO
    ENDIF

    ! set the horizontal attribute pointer of the 3D p_patch
    IF (PRESENT(patch_3d)) THEN
      ALLOCATE(patch_3d, stat=error_status)
      IF (error_status/=success) THEN
        CALL finish(TRIM(routine), 'allocation of patch_3D failed')
      ENDIF
      patch_3d%p_patch_2d => p_patch
    END IF

    CALL message("build_decomposition", "Done.")
  END SUBROUTINE build_decomposition
  !----------------------------------------------------------------------------


  ! -----------------------------------------------------------------
  !
  ! This subroutine computes the child/edge indices of a given patch
  ! "p_p" from the corresponding information on parent indices stored
  ! in the child patch "p_c". This involves parallel comm., since the
  ! child cells and their parent cells may "live" on different PEs.
  !
  ! We also compute the pc_idx here, i.e. the relative ordering of the
  ! child cells/edges in the parent cell/edge. 
  !
  !          *--------------o--------------*
  !           \__   4     _/ \_    2    __/
  !              \__    _/  3  \_    __/
  !                 \_ /_________\ _/
  !                   o           o  
  !                    \__  1  __/
  !                       \_ _/
  !                         *
  !
  ! input:   n_patch_cells, n_patch_cells_g  (on child and parent patch)
  !          n_patch_edges, n_patch_edges_g  (on child and parent patch)
  !          parent_glb_idx/blk  on child patch            (cells and edges)
  !          glb_index           on child and parent patch (cells and edges)
  !          edge_idx/blk        on child patch
  !          decomp_domain, i.e. cell ownership
  !
  ! output:  p_p%cells/edges%child_idx
  !          p_c%cells/edges%pc_idx (if not "is_local_parent")
  !          p_p%cells/edges%child_id
  !          p_p%refin_c_ctrl
  !
  ! 06/2016 : F. Prill, DWD
  !
  SUBROUTINE set_child_indices(description, p_c, p_p, is_local_parent)
    CHARACTER(LEN=*), INTENT(IN) :: description !< description string (for debugging purposes)
    TYPE(t_patch), TARGET, INTENT(inout) :: p_c ! child patch
    TYPE(t_patch), TARGET, INTENT(inout) :: p_p ! parent patch
    LOGICAL, INTENT(IN) :: is_local_parent
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//":set_child_indices"
    INTEGER,      PARAMETER :: UNDEFINED_VALUE = 99
    LOGICAL,      PARAMETER :: EDGE12_REORDER  = .FALSE.

    INTEGER              :: i, j, jc_c, jb_c, communicator, i1, i2, iidx,                    &
      &                     jc_e, jb_e, ordered(4), iglb, je, refin_c
    INTEGER, ALLOCATABLE :: in_data(:), dst_idx(:), out_data(:,:), out_pc_data(:,:),         &
      &                     out_data34(:,:), out_data_e(:,:), out_count(:,:),                &
      &                     out_data_e3(:,:), out_data_c3(:,:), out_child_id(:,:),           &
      &                     je_loc(:,:), out_count_e3(:,:), out_count_e(:,:),                &
      &                     out_pc_count(:,:), out_id_count(:,:)
    LOGICAL              :: lfound, set_pc_idx
    INTEGER, POINTER     :: parent_idx_c(:,:), parent_blk_c(:,:),                            &
      &                     parent_idx_e(:,:), parent_blk_e(:,:)


    IF (msg_level >= 5)  CALL message(routine, TRIM(description)//" - Enter.")

    communicator = p_comm_work

    parent_idx_c => p_c%cells%parent_glb_idx
    parent_blk_c => p_c%cells%parent_glb_blk
    parent_idx_e => p_c%edges%parent_glb_idx
    parent_blk_e => p_c%edges%parent_glb_blk

    set_pc_idx = .NOT. is_local_parent

    ! -----------------------------------------------------------------
    ! --- create CELL indices -----------------------------------------
    ! -----------------------------------------------------------------

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- send the four child cells to the parent cell ----------------

    ALLOCATE(dst_idx(p_c%n_patch_cells), in_data(p_c%n_patch_cells))

    ! get the parent cell indices:
    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      ! loop only over inner domain
      IF (p_c%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE

      iidx = iidx + 1
      dst_idx(iidx) = idx_1d(parent_idx_c(jc_c,jb_c), parent_blk_c(jc_c,jb_c))
      in_data(iidx) = p_c%cells%decomp_info%glb_index(j)
    END DO

    ! consistency check
    IF (MAXVAL(in_data(1:iidx)) > p_c%n_patch_cells_g) THEN
      CALL finish(routine, "Error! - MAXVAL(in_data(1:iidx)) > p_c%n_patch_cells_g")
    END IF

    ALLOCATE(out_data(4,p_p%n_patch_cells), out_count(4,p_p%n_patch_cells))
    out_data(:,:)  = 0
    out_count(:,:) = 0

    ! communicate child index between processors:
    CALL reshuffle("send cell child indices to parent", dst_idx(1:iidx), in_data(1:iidx), &
      &            p_p%cells%decomp_info%glb_index, p_p%n_patch_cells_g, communicator, out_data, out_count)
    DEALLOCATE(in_data, dst_idx, out_count)      

    ! consistency check
    IF (MAXVAL(out_data) > p_c%n_patch_cells_g) THEN
      CALL finish(routine, "Error! - MAXVAL(out_data) > p_c%n_patch_cells_g")
    END IF

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- each child cell sends its three edge indices to the parent
    !     cell. Edges which arrive twice are those of the "central"
    !     child.

    ALLOCATE(dst_idx(3*p_c%n_patch_cells), in_data(3*p_c%n_patch_cells))

    ! get the cell edge indices:
    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      ! loop only over inner domain
      IF (p_c%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE
      DO i=1,3
        jc_e = p_c%cells%edge_idx(jc_c,jb_c,i)
        jb_e = p_c%cells%edge_blk(jc_c,jb_c,i)
        iidx = iidx + 1
        dst_idx(iidx) = idx_1d(parent_idx_c(jc_c,jb_c), parent_blk_c(jc_c,jb_c))
        in_data(iidx) = p_c%edges%decomp_info%glb_index(idx_1d(jc_e,jb_e))
      END DO
    END DO

    ALLOCATE(out_data_e(9,p_p%n_patch_cells), & ! each refined cell has 9 sub-edges
      &      out_count_e(9,p_p%n_patch_cells))
    out_data_e(:,:) = 0
    out_count_e(:,:)  = 0

    ! communicate between processors:
    CALL reshuffle("send cell-edge indices to parent", dst_idx(1:iidx), in_data(1:iidx), &
      &            p_p%cells%decomp_info%glb_index, p_p%n_patch_cells_g, communicator,   &
      &            out_data_e, out_count_e)

    DO j = 1, p_p%n_patch_cells
      IF (SUM(out_count_e(:,j)) == 0)  CYCLE

      ! consistency check: there must be exactly 3 inner edges
      IF (COUNT(out_count_e(:,j) > 1) /= 3) THEN
        WRITE (0,*) "out_count(:,j) = ", out_count_e(:,j)
        WRITE (0,*) "out_data_e(:,j) = ", out_data_e(:,j)
        CALL finish(routine, "edge counting went wrong")
      END IF
    END DO

    DEALLOCATE(in_data, dst_idx)      

    ! --- PARENT PATCH -> CHILD PATCH
    ! --- now, send back the "counting" of the edges to the child patch

    ALLOCATE(dst_idx(12*p_p%n_patch_cells), in_data(12*p_p%n_patch_cells))
    dst_idx = 0
    in_data = 0

    iidx = 0
    DO j = 1, p_p%n_patch_cells
      IF (SUM(out_count_e(:,j)) == 0)  CYCLE
      DO i=1,SIZE(out_data_e,1)
        IF (out_count_e(i,j) > 1) THEN
          DO i1=1,4
            IF (out_data(i1,j) <= 0)  CYCLE
            iidx = iidx + 1
            dst_idx(iidx) = out_data(i1,j)
            in_data(iidx) = out_data_e(i,j)
          END DO
        END IF
      END DO
    END DO

    ! consistency check
    IF (MAXVAL(dst_idx(1:iidx)) > p_c%n_patch_cells_g) THEN
      CALL finish(routine, "Error! - MAXVAL(dst_idx(1:iidx)) > p_c%n_patch_cells_g")
    END IF

    ALLOCATE(out_data_e3(3,p_c%n_patch_cells), out_count_e3(3,p_c%n_patch_cells))
    out_data_e3(:,:)  = 0
    out_count_e3(:,:) = 0

    ! communicate between processors:
    CALL reshuffle("send inner child edges", dst_idx(1:iidx), in_data(1:iidx),         &
      &            p_c%cells%decomp_info%glb_index, p_c%n_patch_cells_g, communicator, &
      &            out_data_e3, out_count_e3)

    DEALLOCATE(in_data, dst_idx,out_data_e, out_count_e, out_count_e3)      

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- determine the global index of the "inner child cell" no. 3
    !     and send this to parent cell:

    ALLOCATE(dst_idx(p_c%n_patch_cells), in_data(p_c%n_patch_cells))
    dst_idx = 0
    in_data = 0

    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      ! loop only over inner domain
      IF (p_c%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE

      lfound = .TRUE.
      EDGE_LOOP : DO i=1,3  ! loop over edges
        jc_e = p_c%cells%edge_idx(jc_c,jb_c,i)
        jb_e = p_c%cells%edge_blk(jc_c,jb_c,i)
        iglb = p_c%edges%decomp_info%glb_index(idx_1d(jc_e,jb_e))
        IF (.NOT. ANY(out_data_e3(:,idx_1d(jc_c,jb_c)) == iglb)) THEN
          lfound = .FALSE. ; EXIT EDGE_LOOP
        END IF
      END DO EDGE_LOOP

      IF (lfound) THEN
        iidx = iidx + 1
        dst_idx(iidx) = idx_1d(parent_idx_c(jc_c,jb_c), parent_blk_c(jc_c,jb_c))
        in_data(iidx) = p_c%cells%decomp_info%glb_index(j)
      END IF
    END DO

    ALLOCATE(out_data_c3(1,p_p%n_patch_cells), out_count(1,p_p%n_patch_cells))
    out_data_c3(:,:) = 0

    ! communicate between processors:
    CALL reshuffle("send child 3 to parent cell", dst_idx(1:iidx), in_data(1:iidx),    &
      &            p_p%cells%decomp_info%glb_index, p_p%n_patch_cells_g, communicator, &
      &            out_data_c3, out_count)

    DEALLOCATE(in_data, dst_idx,out_data_e3, out_count)      

    ! Now, we order the (global) child indices, s.t. 
    ! - the child 3 is in correct position, and
    ! - the remaining children are ordered according to their index.
    DO j = 1, p_p%n_patch_cells
      ! find position of child no. 3
      i1 = 0
      DO i2=1,4
        IF (out_data(i2,j) == out_data_c3(1,j))  i1=i2
      END DO
      ! reordering:
      IF (i1 /= 3)  CALL swap_int(out_data(i1,j), out_data(3,j))
      IF (out_data(1,j) > out_data(2,j))  CALL swap_int(out_data(1,j), out_data(2,j))
      IF (out_data(2,j) > out_data(4,j))  CALL swap_int(out_data(2,j), out_data(4,j))
      IF (out_data(1,j) > out_data(2,j))  CALL swap_int(out_data(1,j), out_data(2,j))
    END DO

    DEALLOCATE(out_data_c3)

    ! --- PARENT PATCH -> CHILD PATCH
    ! --- create CELL pc_idx array ------------------------------------
    !
    ! The "pc_idx" contains the relative ordering of the child cells
    ! in the parent cell.

    ALLOCATE(dst_idx(4*p_p%n_patch_cells), in_data(4*p_p%n_patch_cells))

    iidx = 0
    DO j = 1, p_p%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      ! loop only over inner domain
      IF (p_p%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE
      IF (out_data(1,j) <= 0)  CYCLE

      DO i=1,4
        iidx = iidx + 1
        dst_idx(iidx) = out_data(i,j)
        in_data(iidx) = i
      END DO
    END DO

    ALLOCATE(out_pc_data(1,p_c%n_patch_cells), out_count(1,p_c%n_patch_cells))
    out_pc_data(:,:) = 0

    ! communicate child index between processors:
    CALL reshuffle("send pc_idx to child", dst_idx(1:iidx), in_data(1:iidx), p_c%cells%decomp_info%glb_index, &
      &            p_c%n_patch_cells_g, communicator, out_pc_data, out_count)

    ! --- as a result, we now have
    !     out_pc_data(1,j) = p_c%cells%pc_idx(j)
    !     out_data(i,j) = cells%child_index(j,i)
    IF (set_pc_idx) THEN
      DO j = 1, p_c%n_patch_cells
        jc_c = idx_no(j)
        jb_c = blk_no(j)
        IF (out_pc_data(1,j) <= 0)  CYCLE
        p_c%cells%pc_idx(jc_c,jb_c) = out_pc_data(1,j)
      END DO
    END IF
    DO j = 1, p_p%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      IF (out_data(1,j) <= 0)  CYCLE
      p_p%cells%child_idx(jc_c,jb_c,:) = idx_no( out_data(:,j) )
      p_p%cells%child_blk(jc_c,jb_c,:) = blk_no( out_data(:,j) )
    END DO

    DEALLOCATE(out_pc_data,out_data, in_data, dst_idx, out_count)      

    ! -----------------------------------------------------------------
    ! --- create EDGE indices -----------------------------------------
    ! -----------------------------------------------------------------
    !
    ! The situation is a bit more complicated here: For a given edge,
    ! children 1 and 2 are the sub-triangle edges coinciding with the
    ! parent edge; children 3 and 4 are the edges of the inner
    ! sub-triangles having (roughly) the same orientation as the
    ! parent edge.

    ! --- CHILD PATCH -> PARENT PATCH
    ! --- in a first step, we communicate child edges 3 and 4:
    ALLOCATE(dst_idx(3*p_c%n_patch_cells), in_data(3*p_c%n_patch_cells))

    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      IF (p_c%cells%pc_idx(jc_c,jb_c) == 3) THEN
        DO i=1,3
          jc_e = p_c%cells%edge_idx(jc_c,jb_c,i)
          jb_e = p_c%cells%edge_blk(jc_c,jb_c,i)
          iidx = iidx + 1
          dst_idx(iidx) = idx_1d(parent_idx_e(jc_e,jb_e), parent_blk_e(jc_e,jb_e))
          in_data(iidx) = p_c%edges%decomp_info%glb_index(idx_1d(jc_e,jb_e))
        END DO
      END IF
    END DO

    ALLOCATE(out_data34(2,p_p%n_patch_edges), out_count(2,p_p%n_patch_edges))
    out_data34(:,:) = 0
    out_count(:,:)  = 0

    ! communicate data between processors:
    CALL reshuffle("send child edges 3,4", dst_idx(1:iidx), in_data(1:iidx), p_p%edges%decomp_info%glb_index, &
      &            p_p%n_patch_edges_g, communicator, out_data34, out_count)
    DEALLOCATE(in_data, dst_idx, out_count)      


    ! --- CHILD PATCH -> PARENT PATCH
    ! --- in a second step, we communicate *all four* edge children

    ALLOCATE(dst_idx(p_c%n_patch_edges), in_data(p_c%n_patch_edges))
    iidx = 0
    DO j = 1, p_c%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)
      IF (p_c%edges%decomp_info%decomp_domain(jc_e,jb_e) > 1) CYCLE

      iidx = iidx + 1
      dst_idx(iidx) = idx_1d(parent_idx_e(jc_e,jb_e), parent_blk_e(jc_e,jb_e))
      in_data(iidx) = p_c%edges%decomp_info%glb_index(j)
    END DO

    ALLOCATE(out_data(4,p_p%n_patch_edges), out_count(4,p_p%n_patch_edges))
    out_data(:,:)  = 0
    out_count(:,:) = 0

    ! communicate data between processors:
    CALL reshuffle("send edge indices to parent", dst_idx(1:iidx), in_data(1:iidx), &
      &            p_p%edges%decomp_info%glb_index, p_p%n_patch_edges_g, communicator, out_data, out_count)

    DEALLOCATE(dst_idx, in_data, out_count)

    IF (MAXVAL(out_data) >  p_c%n_patch_edges_g) THEN
      WRITE (0,*) "MAXVAL(out_data)=",MAXVAL(out_data)," > p_c%n_patch_edges_g=",p_c%n_patch_edges_g
      CALL finish(routine, "MAXVAL error 1!")
    END IF

    ! -----------------------------------------------------------------
    ! - build a temporary array which denotes for all parent grid
    !   edges, if the local edge numbering wrt. to the left/right cell
    !   is equal or different.

    IF (EDGE12_REORDER) THEN
      ALLOCATE(dst_idx(3*p_p%n_patch_cells), in_data(3*p_p%n_patch_cells))
      iidx = 0
      DO j = 1, p_p%n_patch_cells
        jc_c = idx_no(j)
        jb_c = blk_no(j)
        DO i=1,3
          ! each edge is traversed twice: 
          jc_e = p_p%cells%edge_idx(jc_c,jb_c,i)
          jb_e = p_p%cells%edge_blk(jc_c,jb_c,i)        
          IF (p_p%edges%decomp_info%decomp_domain(jc_e,jb_e) > 1)  CYCLE
          je = p_p%edges%decomp_info%glb_index(idx_1d(jc_e,jb_e))

          iidx = iidx + 1
          dst_idx(iidx) = je
          in_data(iidx) = i
        END DO
      END DO

      ALLOCATE(je_loc(2,p_p%n_patch_edges), out_count(2,p_p%n_patch_edges))
      je_loc(:,:)    = 0
      out_count(:,:) = 0

      CALL reshuffle("send local edge index", dst_idx(1:iidx), in_data(1:iidx), p_p%edges%decomp_info%glb_index, &
        &            p_p%n_patch_edges_g, communicator, je_loc, out_count)
    END IF

    ! -----------------------------------------------------------------
    ! - copy results

    DO j = 1, p_p%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)

      ! we need to reorder the edge children 3,4 to the end
      ordered(:) = 0
      iidx = 0
      DO i=1,4 
        IF (.NOT. ANY(out_data34(:,j) == out_data(i,j))) THEN
          iidx = iidx + 1
          ordered(iidx) = out_data(i,j)
        END IF
      END DO
      ordered(3:4)  = out_data34(:,j)

      ! -----------------------------------------------------------------
      ! Do some reordering which would not be strictly necessary but
      ! which maintains backward compatibility with the "old" grid
      ! generator.
      ! 
      ! Note that we cannot resemble the numbering in the nest
      ! boundary region here, where the reordering leads to cell and
      ! edge indices which break, e.g., the rule that four consecutive
      ! cell indices form the sub-cells of one parent grid cell.

      ! ascending order of edge index 3,4:
      IF ((ordered(3) > ordered(4)) .AND. (ordered(4) > 0)) THEN
        CALL swap_int(ordered(3), ordered(4))
      END IF

      IF (EDGE12_REORDER) THEN
        !   The rule for edge indices 1,2 is as follows: determine the
        !   local edge index (1/2/3) of the left and right parent grid
        !   cell. If these local edge indices match, then
        !   edges%child_idx(1:2) are ordered in an descending way,
        !   otherwise ascending ordering.
        !
        !   The edge indices 3,4 are always in ascending order.
        !
        ! Behind these rules are three observations:
        !
        ! 1. for a given parent grid edge the left/right cell with the
        !    smaller index defines the child edge numbers, but the
        !    left/right cell with the *larger* index sets the entries
        !    1,2 in the edges%child_idx array (or, more precisely,
        !    overwrites previous settings).
        !
        ! 2. If the local edge index for the parent grid is 1 or 2, then
        !    the indices of the sub-edges are ordered ascendingly, but
        !    if the local edge index for the parent grid is 3, then the
        !    indices are ordered descendingly, (since the child cells
        !    are then #4 and #1).
        !
        ! 3. The edges%child_idx(3) is always set by the *smaller* cell,
        !    and the edges%child_idx(4) by the larger index cell.
        
        ! (do not distinguish between 1/2:)
        !      IF (je_loc(1,j) == 2)  je_loc(1,j) = 1
        IF (je_loc(2,j) == 2)  je_loc(2,j) = 1

        IF (je_loc(2,j) == 0) THEN
          IF (ordered(4) /= 0) THEN
            ! (je_loc(2,j) == 0 can mean that both je_loc values were
            ! equal: descending order)
            IF (ordered(1) < ordered(2))  THEN
              ordered(:) = ordered((/2,1,3,4/))
            END IF
          ELSE
            ! (je_loc(2,j) == 0 can also mean that we are at the
            !  boundary: ascending order
            IF (ordered(1) > ordered(2)) THEN
              ordered(:) = ordered((/2,1,3,4/))
            END IF
          END IF
        ELSE
          ! standard case: ascending order, if je_loc not matching,
          ! (the other option - ascending ordering for matching - was handled above):
          IF ((je_loc(1,j) /= je_loc(2,j)) .AND. (ordered(1) > ordered(2))) THEN
            ordered(:) = ordered((/2,1,3,4/))
          END IF
        END IF
      ELSE
        IF (ordered(1) > ordered(2))  THEN ! DEVELOPMENT: Simple ascending ordering
          ordered(:) = ordered((/2,1,3,4/))
        END IF
      END IF

      out_data(:,j) = ordered(:)
    END DO
    IF (EDGE12_REORDER) THEN
      DEALLOCATE(dst_idx, in_data, out_data34, je_loc, out_count)
    END IF


    ! --- PARENT PATCH -> CHILD PATCH
    ! --- create EDGE pc_idx array ------------------------------------
    !
    ! The "pc_idx" contains the relative ordering of the child edges (1...4).

    ALLOCATE(dst_idx(4*p_p%n_patch_edges), in_data(4*p_p%n_patch_edges))

    iidx = 0
    DO j = 1, p_p%n_patch_edges
      DO i=1,4
        IF (out_data(i,j) > 0) THEN
          iidx = iidx + 1
          dst_idx(iidx) = out_data(i,j)
          in_data(iidx) = i
        END IF
      END DO
    END DO

    ALLOCATE(out_pc_data(1,p_c%n_patch_edges), out_pc_count(1,p_c%n_patch_edges))
    out_pc_data(:,:)  = 0
    out_pc_count(:,:) = 0

    ! communicate child index between processors:
    CALL reshuffle("send edge pc_idx", dst_idx(1:iidx), in_data(1:iidx), p_c%edges%decomp_info%glb_index, &
      &            p_c%n_patch_edges_g, communicator, out_pc_data, out_pc_count)

    ! --- as a result, we now have
    !     out_pc_data(1,j) = p_c%edges%pc_idx(j)
    !       and 
    !     out_data(:,j) = p_p%edges%child_idx(j,:)
    IF (set_pc_idx) THEN
      DO j = 1, p_c%n_patch_edges
        jc_e = idx_no(j)
        jb_e = blk_no(j)
        p_c%edges%pc_idx(jc_e,jb_e) = out_pc_data(1,j)
      END DO
    END IF
    DO j = 1, p_p%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)
      IF (ALL(out_data(:,j) == 0))  CYCLE
      p_p%edges%child_idx(jc_e,jb_e,:) = idx_no( out_data(:,j) )
      p_p%edges%child_blk(jc_e,jb_e,:) = blk_no( out_data(:,j) )
    END DO

    DEALLOCATE(out_data, in_data, dst_idx, out_pc_count)      

    ! -----------------------------------------------------------------
    ! --- create CHILD ID's -------------------------------------------
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    ! --- calculate cells%child_id
    !
    ! To calculate the child_id for each child cell, we loop over
    ! the cells of the child patch and collect the global indices of
    ! their parent cells. Then we send the ID of the child domain to
    ! these parent cells. 

    ALLOCATE(dst_idx(p_c%n_patch_cells), in_data(p_c%n_patch_cells))

    ! get the parent cell indices:
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)
      dst_idx(j) = idx_1d(parent_idx_c(jc_c,jb_c), parent_blk_c(jc_c,jb_c))
    END DO
    ! set the cell child_id for the (global) parent grid
    in_data(:) = p_c%id

    ALLOCATE(out_child_id(1,p_p%n_patch_cells), out_id_count(1,p_p%n_patch_cells))
    out_child_id(:,:) = -1
    out_id_count(:,:) = -1

    ! communicate ID data between processors:
    CALL reshuffle("send cell id to parent", dst_idx, in_data, p_p%cells%decomp_info%glb_index, &
      &            p_p%n_patch_cells_g, communicator, out_child_id, out_id_count)

    ! communication finished. now copy the result to the local arrays:
    DO j = 1, p_p%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      IF (out_child_id(1,j) < 0)  CYCLE

      ! check for overlap (i.e.: child id already exists)
      IF (p_p%cells%child_id(jc_c,jb_c) == 0) THEN
        p_p%cells%child_id(jc_c,jb_c) = out_child_id(1,j)
      ELSE
        IF (p_p%cells%child_id(jc_c,jb_c) /= out_child_id(1,j)) THEN
          CALL finish(routine, "Cell in DOM "//int2string(p_c%id,'(i0)')//" overlaps with another domain!")
        END IF
      END IF
    END DO
    DEALLOCATE(out_child_id, in_data, dst_idx, out_id_count)

    ! -----------------------------------------------------------------
    ! --- calculate edges%child_id
    !
    ! To calculate the child_id for each child edge, we loop over the
    ! edges of the child domain and collect the global indices of
    ! their parent edges. Then we send the ID of the child domain to
    ! these parent edges. 

    ALLOCATE(dst_idx(p_c%n_patch_edges), in_data(p_c%n_patch_edges))

    ! get the parent cell indices:
    DO j = 1, p_c%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)
      dst_idx(j) = idx_1d(parent_idx_e(jc_e,jb_e), parent_blk_e(jc_e,jb_e))
    END DO
    ! set the edge child_id for the (global) parent grid:
    in_data(:) = p_c%id

    ALLOCATE(out_child_id(1,p_p%n_patch_edges), out_id_count(1,p_p%n_patch_edges))
    out_child_id(:,:) = -1
    out_id_count(:,:)    = -1

    ! communicate ID data between processors:
    CALL reshuffle("send edge ID to parent", dst_idx, in_data, p_p%edges%decomp_info%glb_index, &
      &            p_p%n_patch_edges_g, communicator, out_child_id, out_id_count)

    ! communication finished. now copy the result to the local arrays:
    DO j = 1, p_p%n_patch_edges
      jc_e = idx_no(j)
      jb_e = blk_no(j)

      IF (out_child_id(1,j) < 0)  CYCLE

      ! check for overlap (i.e.: child id already exists)
      IF (p_p%edges%child_id(jc_e,jb_e) == 0) THEN
        p_p%edges%child_id(jc_e,jb_e) = out_child_id(1,j)
      ELSE
        IF (p_p%edges%child_id(jc_e,jb_e) /= out_child_id(1,j)) THEN
          CALL finish(routine, "Edge in DOM "//int2string(p_c%id,'(i0)')//" overlaps with another domain!")
        END IF
      END IF
    END DO
    DEALLOCATE(out_child_id, in_data, dst_idx, out_id_count)


    ! -----------------------------------------------------------------
    ! --- create REFIN_C_CTRL of parent grid --------------------------
    ! -----------------------------------------------------------------

    ALLOCATE(dst_idx(p_c%n_patch_cells), in_data(p_c%n_patch_cells))
  
    ! send the children's refin_ctrl to the parent cell:
    iidx = 0
    DO j = 1, p_c%n_patch_cells
      jc_c = idx_no(j) ; jb_c = blk_no(j)
      ! loop only over inner domain
      IF (p_c%cells%decomp_info%decomp_domain(jc_c,jb_c) /= 0)  CYCLE
      IF (p_c%cells%refin_ctrl(jc_c,jb_c) <= 0)  CYCLE

      refin_c = (p_c%cells%refin_ctrl(jc_c,jb_c)+1)/2 

      iidx = iidx + 1
      dst_idx(iidx) = idx_1d(parent_idx_c(jc_c,jb_c), parent_blk_c(jc_c,jb_c))
      in_data(iidx) = refin_c
    END DO
    
    ALLOCATE(out_data(4,p_p%n_patch_cells), out_count(4,p_p%n_patch_cells))
    out_data(:,:)  = 0
    out_count(:,:) = 0
 
    ! communicate between processors:
    CALL reshuffle("send refin_c_ctrl to parent", dst_idx(1:iidx), in_data(1:iidx),    &
      &            p_p%cells%decomp_info%glb_index, p_p%n_patch_cells_g, communicator, &
      &            out_data, out_count)
  
    ! communication finished. now copy the result to the local arrays:
    DO j = 1, p_p%n_patch_cells
      jc_c = idx_no(j)
      jb_c = blk_no(j)

      ! first, set all parent cell refin_ctrl flags to min_rlcell_int
      ! based on the information in "p_p%cells%child_id": this handles
      ! the innermost part of the domain:
      IF (p_p%cells%child_id(jc_c,jb_c) == p_c%id) THEN
        refin_c = min_rlcell_int
      ELSE
        refin_c = UNDEFINED_VALUE
      END IF

      IF (out_data(1,j) /= 0) THEN
        refin_c = MAX( -1*out_data(1,j), min_rlcell_int )
      END IF

      IF (refin_c /= UNDEFINED_VALUE) THEN
        p_p%cells%refin_ctrl(jc_c,jb_c) = refin_c
      END IF
    END DO
    DEALLOCATE(out_data, in_data, dst_idx, out_count)

    IF (msg_level >= 5)  CALL message(routine, TRIM(description)//" - Done.")

  CONTAINS

    SUBROUTINE swap_int(a,b)
      INTEGER,  INTENT(INOUT) :: a,b
      ! local variables
      INTEGER :: tmp
      tmp = a ; a = b; b = tmp
    END SUBROUTINE swap_int

  END SUBROUTINE set_child_indices

 
END MODULE mo_build_decomposition
