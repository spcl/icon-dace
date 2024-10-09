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

! Basic tools for the patch

!----------------------------
#include "omp_definitions.inc"
!----------------------------

MODULE mo_grid_tools
  !-------------------------------------------------------------------------
  USE mo_kind,               ONLY: wp
  USE mo_exception,          ONLY: finish, warning
  USE mo_parallel_config,    ONLY: nproma
  USE mo_model_domain,       ONLY: t_patch, t_patch_3D
  USE mo_math_utilities,     ONLY: gvec2cvec, gc2cc
  USE mo_math_types
  USE mo_grid_subset,        ONLY: t_subset_range, get_index_range, t_subset_indexed, &
    & block_no, index_no, index_1d
  USE mo_impl_constants,     ONLY: on_edges, land
  USE mo_mpi,                ONLY: get_my_mpi_work_id, p_sum, get_my_mpi_work_communicator, &
    & my_process_is_stdio, work_mpi_barrier

  USE mo_decomposition_tools,ONLY: t_glb2loc_index_lookup, get_local_index, get_valid_local_index
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC :: calculate_patch_cartesian_positions
  PUBLIC :: rescale_grid
  PUBLIC :: calculate_edge_area
  PUBLIC :: get_oriented_edges_from_global_vertices
  PUBLIC :: find_oriented_edge_from_vertices
  PUBLIC :: create_dummy_cell_closure
  PUBLIC :: check_global_indexes
  
  
  !-------------------------------------------------------------------------
  
CONTAINS

  !-------------------------------------------------------------------------
  !>
  ! Fills the edges subset defined by the global_vertex_array.
  ! Also returns the relative orientation
   SUBROUTINE get_oriented_edges_from_global_vertices(edge_subset, orientation, patch_2D, patch_3D, &
     & global_vertex_array, subset_name)
    TYPE(t_subset_indexed), INTENT(inout) :: edge_subset
    REAL(wp), POINTER :: orientation(:)
    TYPE(t_patch), TARGET, INTENT(in),    OPTIONAL :: patch_2D  ! nag does not return the values in subset
    TYPE(t_patch_3D), TARGET, INTENT(in), OPTIONAL :: patch_3D  ! nag does not return the values in subset
    INTEGER :: global_vertex_array(:)   ! intent in
    CHARACTER(len=*), OPTIONAL :: subset_name

    INTEGER :: my_proc_id, i, max_allocation_size, owned_edges
    INTEGER :: edge_index, edge_block
    INTEGER :: local_vertex_1d_index(2), vertex_index(2), vertex_block(2)
    INTEGER, ALLOCATABLE :: tmp_edge_block_array(:), tmp_edge_index_array(:)
    REAL(wp) :: edge_orientation
    REAL(wp), ALLOCATABLE :: tmp_orientation(:)
    INTEGER, POINTER :: local_vertex_array(:), owner_edge_local(:)
    TYPE(t_glb2loc_index_lookup), POINTER :: verts_glb2loc_index

    INTEGER :: total_subset_edges
    INTEGER :: my_mpi_work_communicator
    CHARACTER(*), PARAMETER :: method_name = "mo_grid_tools:get_oriented_edges_from_global_vertices"

    edge_subset%size               = 0
    edge_subset%recommended_stride = 0
    edge_subset%entity_location    = on_edges
    NULLIFY(edge_subset%patch)
    NULLIFY(edge_subset%patch_3D)
    ! NULLIFY(edge_subset%block)
    ! NULLIFY(edge_subset%idx)
    NULLIFY(orientation)
    edge_subset%name=""
    IF (PRESENT(subset_name)) edge_subset%name=subset_name

    IF (PRESENT(patch_3D)) THEN
      edge_subset%patch_3D => patch_3D
      edge_subset%patch    => patch_3D%p_patch_2D(1)
    ELSEIF (PRESENT(patch_2D)) THEN
      edge_subset%patch    => patch_2D
    ELSE
      CALL finish(method_name, "patch is not present")
    ENDIF
    verts_glb2loc_index => edge_subset%patch%verts%decomp_info%glb2loc_index
    owner_edge_local    => edge_subset%patch%edges%decomp_info%owner_local

    my_proc_id = get_my_mpi_work_id()

    ! temporary array for keeping track of what's local
    max_allocation_size = SIZE(global_vertex_array)
    DO i=1, max_allocation_size
      IF ( global_vertex_array(i) <= 0 ) THEN
        max_allocation_size = i - 1
        EXIT
      ENDIF
    ENDDO

    ALLOCATE(tmp_edge_block_array(max_allocation_size), &
      & tmp_edge_index_array(max_allocation_size),      &
      & tmp_orientation(max_allocation_size))
    owned_edges = 0
    IF (max_allocation_size > 1) &
      & local_vertex_1d_index(1) = get_local_index(verts_glb2loc_index, &
      &                                            global_vertex_array(1))

!     write(0,*) 1, ":", global_vertex_array(1),  local_vertex_1d_index(1)
    DO i=2, max_allocation_size
      ! get a pair of vertices
      local_vertex_1d_index(2) = get_local_index(verts_glb2loc_index, &
        &                                        global_vertex_array(i))

      IF (local_vertex_1d_index(1) > 0 .AND. local_vertex_1d_index(2) > 0 ) THEN
        ! find the edge and orientation of these vertices
        vertex_block(1) = block_no(local_vertex_1d_index(1))
        vertex_index(1) = index_no(local_vertex_1d_index(1))
        vertex_block(2) = block_no(local_vertex_1d_index(2))
        vertex_index(2) = index_no(local_vertex_1d_index(2))

        CALL find_oriented_edge_from_vertices(patch=edge_subset%patch, &   ! in
          & vertex_blocks=vertex_block, vertex_indexes=vertex_index,   &   ! in
          & edge_block=edge_block, edge_index=edge_index,              &   ! out
          & edge_orientation=edge_orientation)                             ! out

        IF (edge_index > 0) THEN
!           write(0,*) "Found edge :", edge_block, edge_index
          IF (owner_edge_local(index_1d(idx=edge_index, block=edge_block)) == my_proc_id) THEN

            write(0,*) i-1, ":", global_vertex_array(i-1),  local_vertex_1d_index(1), &
              & edge_subset%patch%verts%decomp_info%glb_index(local_vertex_1d_index(1))
            write(0,*) i, ":", global_vertex_array(i),  local_vertex_1d_index(2), &
              & edge_subset%patch%verts%decomp_info%glb_index(local_vertex_1d_index(2))
          
            write(0,*) my_proc_id, " owns ", global_vertex_array(i-1), global_vertex_array(i), &
             & edge_subset%patch%edges%decomp_info%glb_index((edge_block-1)*nproma+edge_index), &
             & " halo_level=", edge_subset%patch%edges%decomp_info%halo_level(edge_index, edge_block)
             
            owned_edges = owned_edges + 1
            tmp_edge_block_array(owned_edges) = edge_block
            tmp_edge_index_array(owned_edges) = edge_index
            tmp_orientation(owned_edges)      = edge_orientation
          ENDIF
        ENDIF
      ENDIF

      ! get next pair
      local_vertex_1d_index(1) = local_vertex_1d_index(2)
    ENDDO

    !now fill the edge_subset if not empty
    IF (owned_edges > 0) THEN

      edge_subset%size = owned_edges
      edge_subset%recommended_stride = 1  ! needs to be calculated
      ALLOCATE(edge_subset%block(owned_edges), &
        &      edge_subset%idx(owned_edges),   &
        &      orientation(owned_edges))

      DO i=1, owned_edges
        edge_subset%block(i) = tmp_edge_block_array(i)
        edge_subset%idx(i)   = tmp_edge_index_array(i)
        orientation(i)       = tmp_orientation(i)
      ENDDO

    ENDIF

    DEALLOCATE(tmp_edge_block_array, tmp_edge_index_array, tmp_orientation)

    !check the total subset size
    CALL work_mpi_barrier()
    my_mpi_work_communicator = get_my_mpi_work_communicator()
    total_subset_edges = p_sum(edge_subset%size, comm=my_mpi_work_communicator)
    IF (my_process_is_stdio()) THEN
       write(0,*) TRIM(edge_subset%name), " Subset total entities=", total_subset_edges, " in vertices=", max_allocation_size
    ENDIF
    
    CALL work_mpi_barrier()
    
  END SUBROUTINE get_oriented_edges_from_global_vertices
  !----------------------------------------------------

  !-------------------------------------------------------------------------
  !> Returns a grid edge -if exists- and the relative orientation defined by two vertices
  ! Returns 0 if an edge was not found
  SUBROUTINE find_oriented_edge_from_vertices(patch, vertex_blocks, vertex_indexes, &
    & edge_block, edge_index, edge_orientation)
    TYPE(t_patch), TARGET, INTENT(in) :: patch
    INTEGER, INTENT(in)   :: vertex_blocks(2), vertex_indexes(2)
    INTEGER, INTENT(out)  :: edge_block, edge_index
    REAL(wp), INTENT(out) :: edge_orientation
  
    INTEGER :: i, j, edge1_block, edge1_idx
    LOGICAL :: found
    CHARACTER(*), PARAMETER :: method_name = "mo_grid_tools:find_oriented_edge_from_vertices"

    edge_block = 0
    edge_index = 0
    edge_orientation = 0.0_wp

    IF (vertex_indexes(1) == vertex_indexes(2) .AND. vertex_blocks(1) == vertex_blocks(2)) THEN
      CALL finish(method_name, "Identical vertices")
    ENDIF

   ! just brute force search
    found = .false.
    DO i=1, patch%verts%num_edges(vertex_indexes(1), vertex_blocks(1))
      edge1_block = patch%verts%edge_blk(vertex_indexes(1), vertex_blocks(1), i)
      edge1_idx   = patch%verts%edge_idx(vertex_indexes(1), vertex_blocks(1), i)

      DO j=1, patch%verts%num_edges(vertex_indexes(2), vertex_blocks(2))

         IF (patch%verts%edge_idx(vertex_indexes(2), vertex_blocks(2), j) == edge1_idx .AND. &
          &  patch%verts%edge_blk(vertex_indexes(2), vertex_blocks(2), j) == edge1_block) THEN
            found = .true.
          EXIT
        ENDIF

      ENDDO

      IF (found) EXIT

    ENDDO

    IF (found) THEN

      edge_block = edge1_block
      edge_index = edge1_idx
      ! fill the orientation
      IF (patch%edges%vertex_idx(edge_index, edge_block, 1) == vertex_indexes(1) .AND. &
        & patch%edges%vertex_blk(edge_index, edge_block, 1) == vertex_blocks(1)  .AND. &
        & patch%edges%vertex_idx(edge_index, edge_block, 2) == vertex_indexes(2) .AND. &
        & patch%edges%vertex_blk(edge_index, edge_block, 2) == vertex_blocks(2)) THEN
        edge_orientation = patch%edges%tangent_orientation(edge_index, edge_block)
      ELSEIF(patch%edges%vertex_idx(edge_index, edge_block, 1) == vertex_indexes(2) .AND. &
        & patch%edges%vertex_blk(edge_index, edge_block, 1) == vertex_blocks(2)  .AND. &
        & patch%edges%vertex_idx(edge_index, edge_block, 2) == vertex_indexes(1) .AND. &
        & patch%edges%vertex_blk(edge_index, edge_block, 2) == vertex_blocks(1)) THEN
        edge_orientation = -patch%edges%tangent_orientation(edge_index, edge_block)
      ELSE
        CALL finish(method_name, "Edge-Vertex connectivity inconsistency")
      ENDIF

!      write(0,*) get_my_mpi_work_id(), ": edge found:", edge_index, edge_block, edge_orientation
!    ELSE
!      CALL finish("find_oriented_edge_from_vertices","edge not found")


    ENDIF

  END SUBROUTINE find_oriented_edge_from_vertices
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Rescale grids
  ! Note: this does not rescale the cartesian coordinates (used in torus)
  SUBROUTINE rescale_grid( patch, grid_length_rescale_factor)
    
    TYPE(t_patch), INTENT(inout), TARGET ::  patch  ! patch data structure
    REAL(wp), INTENT(in) :: grid_length_rescale_factor
    
    REAL(wp) :: grid_area_rescale_factor
    !-----------------------------------------------------------------------
    grid_area_rescale_factor   = grid_length_rescale_factor * grid_length_rescale_factor
    
    patch%cells%area(:,:)               = patch%cells%area(:,:)               * grid_area_rescale_factor
    patch%verts%dual_area(:,:)          = patch%verts%dual_area(:,:)          * grid_area_rescale_factor
    patch%edges%primal_edge_length(:,:) = patch%edges%primal_edge_length(:,:) * grid_length_rescale_factor
    patch%edges%dual_edge_length(:,:)   = patch%edges%dual_edge_length(:,:)   * grid_length_rescale_factor
    patch%edges%edge_cell_length(:,:,:) = patch%edges%edge_cell_length(:,:,:) * grid_length_rescale_factor
    patch%edges%edge_vert_length(:,:,:) = patch%edges%edge_vert_length(:,:,:) * grid_length_rescale_factor

    ! rescale geometry parameters
    patch%geometry_info%mean_edge_length           = patch%geometry_info%mean_edge_length           * grid_length_rescale_factor
    patch%geometry_info%mean_cell_area             = patch%geometry_info%mean_cell_area             * grid_area_rescale_factor
    patch%geometry_info%mean_dual_edge_length      = patch%geometry_info%mean_dual_edge_length      * grid_length_rescale_factor
    patch%geometry_info%mean_dual_cell_area        = patch%geometry_info%mean_dual_cell_area        * grid_area_rescale_factor
    patch%geometry_info%domain_length              = patch%geometry_info%domain_length              * grid_length_rescale_factor
    patch%geometry_info%domain_height              = patch%geometry_info%domain_height              * grid_length_rescale_factor
    patch%geometry_info%sphere_radius              = patch%geometry_info%sphere_radius              * grid_length_rescale_factor
    patch%geometry_info%mean_characteristic_length = patch%geometry_info%mean_characteristic_length * grid_length_rescale_factor

    IF (patch%geometry_info%mean_characteristic_length == 0.0_wp) &
      & CALL finish("rescale_grid", "mean_characteristic_length=0")
    
  END SUBROUTINE rescale_grid
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! This routine calculates the edge area of the patch
  !!
  SUBROUTINE calculate_edge_area(patch)
    TYPE(t_patch), TARGET, INTENT(inout) :: patch

    INTEGER                       :: je, jb, start_edge_index, end_edge_index
    TYPE(t_subset_range), POINTER :: all_edges
    !-------------------------------------------------------------------------
    all_edges => patch%edges%all

    ! a) the control volume associated to each edge is defined as the
    ! quadrilateral whose edges are the primal edge and the associated dual edge
    !----------------------------------------------------------------------------
    ! loop over all blocks and edges

!ICON_OMP_PARALLEL DO PRIVATE(start_edge_index, end_edge_index,je) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, start_edge_index, end_edge_index)
      DO je = start_edge_index, end_edge_index
        patch%edges%area_edge(je,jb) =  &
            &    patch%edges%primal_edge_length(je,jb)  &
            &  * patch%edges%dual_edge_length(je,jb)
      END DO
    END DO
!ICON_OMP_END_PARALLEL DO

  END SUBROUTINE calculate_edge_area
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  !>
  !! This method_name calculates the cartesian positions stored in patch.
  !!
  !! These values normally are read from the grid files,
  !! this routine is used only for old grids
  !!
  !! Note: cartesian_dual_middle is used only from the ocean,
  !!   which is using only new grids, therefor it is NOT calculated here
  !!
  SUBROUTINE calculate_patch_cartesian_positions ( patch )
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: patch  ! patch on a specific level
    
    TYPE(t_cartesian_coordinates) :: z_vec
    REAL(wp) :: z_lon, z_lat, z_u, z_v  ! location and components of normal
    REAL(wp) :: z_norm                  ! norm of Cartesian normal
    
    INTEGER :: start_index, end_index
    INTEGER :: jb, je                   ! loop indices    
    !-----------------------------------------------------------------------
                
!ICON_OMP_PARALLEL
    ! calculate cells cartesian positions
!ICON_OMP_DO PRIVATE(jb,je, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = patch%cells%all%start_block, patch%cells%all%end_block
      CALL get_index_range(patch%cells%all, jb, start_index, end_index)
      DO je = start_index, end_index            
        ! location of cell center
        patch%cells%cartesian_center(je,jb) = gc2cc(patch%cells%center(je,jb))
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT
    
    ! calculate verts cartesian positions
!ICON_OMP_DO PRIVATE(jb,je, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = patch%verts%all%start_block, patch%verts%all%end_block
      CALL get_index_range(patch%verts%all, jb, start_index, end_index)
      DO je = start_index, end_index            
        ! location of cell center
        patch%verts%cartesian(je,jb) = gc2cc(patch%verts%vertex(je,jb))
      ENDDO
    ENDDO
!ICON_OMP_END_DO NOWAIT

    ! calculate edges positions
!ICON_OMP_DO PRIVATE(jb,je,z_lon,z_lat,z_u,z_v,z_norm,z_vec, start_index, end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = patch%edges%all%start_block, patch%edges%all%end_block
      CALL get_index_range(patch%edges%all, jb, start_index, end_index)
      DO je = start_index, end_index
            
        ! location of edge midpoint
        patch%edges%cartesian_center(je,jb) = gc2cc(patch%edges%center(je,jb))
        z_lon = patch%edges%center(je,jb)%lon
        z_lat = patch%edges%center(je,jb)%lat
        
        ! zonal and meridional component of primal normal
        z_u = patch%edges%primal_normal(je,jb)%v1
        z_v = patch%edges%primal_normal(je,jb)%v2
        
        ! calculate Cartesian components of primal normal
        CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )
        
        ! compute unit normal to edge je
        z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
        z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)
        
        ! save the values in the according type structure of the patch
!         WRITE(*,*) "----------------------------------"
!         write(*,*) "primal_cart_normal:", patch%edges%primal_cart_normal(je,jb)%x(:)
!         write(*,*) "z_vec:", z_vec%x
!         WRITE(*,*) "----------------------------------"
!         IF ( MAXVAL(ABS(patch%edges%primal_cart_normal(je,jb)%x - z_vec%x)) > 0.0001_wp) &
!           CALL finish("","primal_cart_normal(je,jb)%x /=  z_vec%x")
               
        patch%edges%primal_cart_normal(je,jb)%x(1) = z_vec%x(1)
        patch%edges%primal_cart_normal(je,jb)%x(2) = z_vec%x(2)
        patch%edges%primal_cart_normal(je,jb)%x(3) = z_vec%x(3)
        
        ! zonal and meridional component of dual normal
        z_u = patch%edges%dual_normal(je,jb)%v1
        z_v = patch%edges%dual_normal(je,jb)%v2
        
        ! calculate Cartesian components of dual normal
        CALL gvec2cvec( z_u, z_v, z_lon, z_lat, z_vec%x(1), z_vec%x(2), z_vec%x(3) )
        
        ! compute unit normal to edge je
        z_norm = SQRT( DOT_PRODUCT(z_vec%x(1:3),z_vec%x(1:3)) )
        z_vec%x(1:3) = 1._wp / z_norm * z_vec%x(1:3)
        
!         WRITE(*,*) "----------------------------------"
!         write(*,*) "dual_cart_normal:", patch%edges%dual_cart_normal(je,jb)%x(:)
!         write(*,*) "z_vec:", z_vec%x
!         WRITE(*,*) "----------------------------------"
!         IF ( MAXVAL(ABS(patch%edges%dual_cart_normal(je,jb)%x - z_vec%x)) > 0.0001_wp) &
!           CALL finish("","dual_cart_normal(je,jb)%x /=  z_vec%x")
        
        ! save the values in the according type structure of the patch
        patch%edges%dual_cart_normal(je,jb)%x(1) = z_vec%x(1)
        patch%edges%dual_cart_normal(je,jb)%x(2) = z_vec%x(2)
        patch%edges%dual_cart_normal(je,jb)%x(3) = z_vec%x(3)
        
      END DO
      
    END DO
!ICON_OMP_END_DO NOWAIT
!ICON_OMP_END_PARALLEL
    
  END SUBROUTINE calculate_patch_cartesian_positions
  !-------------------------------------------------------------------------

    
  !-------------------------------------------------------------------------
  !>
  ! adds a dummy cell at the end of the existing cells
  ! and creates a "dummy" connectivity for edges and cells that have no neighbor
  ! Note that the dummy cell has to be already allocated!
  ! This is indented for the ocean, with no land cells
  !
  ! The subsets must have been filled in order in order to call this routine.
  SUBROUTINE create_dummy_cell_closure( patch_3D )
    TYPE(t_patch_3D), INTENT(inout), TARGET ::  patch_3D

    TYPE(t_patch), POINTER ::  patch
    INTEGER :: block, idx, start_idx, end_idx, neighbor
    INTEGER :: dummy_cell_blk, dummy_cell_idx

    CHARACTER(LEN=*), PARAMETER :: method_name = 'mo_grid_tools:create_dummy_cell_closure'

    patch => patch_3D%p_patch_2D(1)
    ! write(0,*) "------------ start create_dummy_cell_closure ---------------"
    ! write(0,*) "patch%cells%max_connectivity=", patch%cells%max_connectivity

    dummy_cell_idx = patch%cells%dummy_cell_index
    dummy_cell_blk = patch%cells%dummy_cell_block

   ! Fill edge connectivity with the dummy cell
!ICON_OMP_PARALLEL
!ICON_OMP_DO PRIVATE(block, idx, neighbor, start_idx, end_idx)
    DO block = patch%edges%all%start_block, patch%edges%all%end_block
      CALL get_index_range(patch%edges%all, block, start_idx, end_idx)
      DO idx = start_idx, end_idx
        DO neighbor=1,2
          IF ( patch%edges%cell_idx(idx, block, neighbor) == 0) THEN
            patch%edges%cell_blk(idx, block, neighbor) = dummy_cell_blk
            patch%edges%cell_idx(idx, block, neighbor) = dummy_cell_idx
          ENDIF
        END DO
      END DO
    END DO
!ICON_OMP_END_DO

   ! Fill cell connectivity with the dummy cell
!ICON_OMP_DO PRIVATE(block, idx, neighbor, start_idx, end_idx)
    DO block = patch%cells%all%start_block, patch%cells%all%end_block
      CALL get_index_range(patch%cells%all, block, start_idx, end_idx)
      DO idx = start_idx, end_idx
         IF (block /= dummy_cell_blk .OR. idx /= dummy_cell_idx) THEN ! this is not needed
          ! this is not the dummy cell
          DO neighbor=1, patch%cells%max_connectivity
            IF ( patch%cells%neighbor_idx(idx, block, neighbor) == 0 .OR. &
              & patch%cells%neighbor_blk(idx, block, neighbor) == 0 ) THEN
              patch%cells%neighbor_blk(idx, block, neighbor) = dummy_cell_blk
              patch%cells%neighbor_idx(idx, block, neighbor) = dummy_cell_idx
  !             write(0,*) "Replaced neighbor at ", block, idx, neighbor
            ENDIF
          END DO
         ENDIF
      END DO
    END DO
!ICON_OMP_END_DO
!ICON_OMP_END_PARALLEL

    ! make dummy land
    patch_3D%lsm_c(dummy_cell_idx,:,dummy_cell_blk)                       = land
    patch_3D%surface_cell_sea_land_mask(dummy_cell_idx,dummy_cell_blk)    = land
    patch_3D%wet_c(dummy_cell_idx,:,dummy_cell_blk)                       = 0.0_wp
    patch_3D%wet_halo_zero_c(dummy_cell_idx,:,dummy_cell_blk)             = 0.0_wp
    patch_3D%basin_c(dummy_cell_idx,dummy_cell_blk)                       = 0
    patch_3D%regio_c(dummy_cell_idx,dummy_cell_blk)                       = 0
    patch_3D%bottom_thick_c(dummy_cell_idx,dummy_cell_blk)                = 0.0_wp
    patch_3D%column_thick_c(dummy_cell_idx,dummy_cell_blk)                = 0.0_wp

    patch_3D%p_patch_1D(1)%dolic_c(dummy_cell_idx,dummy_cell_blk)         = 0
    patch_3D%p_patch_1D(1)%prism_thick_c(dummy_cell_idx,:,dummy_cell_blk) = 0.0_wp
    patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(dummy_cell_idx,:,dummy_cell_blk) = 0.0_wp

  END SUBROUTINE create_dummy_cell_closure
  !-----------------------------------------------------------------------
    
  !-------------------------------------------------------------------------
  SUBROUTINE check_global_indexes(patch_2D)
    TYPE(t_patch), TARGET, INTENT(in),    OPTIONAL :: patch_2D  ! nag does not return the values in subset

    INTEGER :: my_proc_id
    INTEGER :: global_index, glb2loc, valid_glb2loc, loc2glb
    INTEGER :: blockNo, indexNo, noOfEdges
    INTEGER :: my_mpi_work_communicator
    TYPE(t_glb2loc_index_lookup), POINTER :: glb2loc_index
    INTEGER, POINTER :: owner_local(:)

    my_proc_id = get_my_mpi_work_id()

    glb2loc_index  => patch_2D%edges%decomp_info%glb2loc_index
    owner_local    => patch_2D%edges%decomp_info%owner_local
!     global_index = 0
!     glb2loc = get_local_index(glb2loc_index, &
!       &  global_index)
!     valid_glb2loc = get_valid_local_index(glb2loc_index, global_index, .TRUE.)
!     write(0,*) my_proc_id, ": edge:", global_index, glb2loc, valid_glb2loc
    
    global_index = 20834
    glb2loc = get_local_index(glb2loc_index, &
      &  global_index)
    valid_glb2loc = get_valid_local_index(glb2loc_index, global_index, .TRUE.)
    IF (glb2loc > 0)  THEN
      loc2glb = patch_2D%edges%decomp_info%glb_index(glb2loc)
      blockNo = block_no(glb2loc)
      indexNo = index_no(glb2loc)
      write(0,*) my_proc_id, ": edge: (", patch_2D%edges%center(indexNo,blockNo)%lon, &
        & patch_2D%edges%center(indexNo, blockNo)%lat, ")", &
        & global_index, glb2loc, valid_glb2loc, loc2glb, owner_local(glb2loc)
    ENDIF

    global_index = 20835
    glb2loc = get_local_index(glb2loc_index, &
      &  global_index)
    valid_glb2loc = get_valid_local_index(glb2loc_index, global_index, .TRUE.)
    IF (glb2loc > 0)  THEN
      loc2glb = patch_2D%edges%decomp_info%glb_index(glb2loc)
      blockNo = block_no(glb2loc)
      indexNo = index_no(glb2loc)
      write(0,*) my_proc_id, ": edge: (", patch_2D%edges%center(indexNo,blockNo)%lon, &
        & patch_2D%edges%center(indexNo, blockNo)%lat, ")", &
        & global_index, glb2loc, valid_glb2loc, loc2glb, owner_local(glb2loc)
    ENDIF
    
    glb2loc_index  => patch_2D%verts%decomp_info%glb2loc_index
    owner_local    => patch_2D%verts%decomp_info%owner_local

    global_index = 7288
    glb2loc = get_local_index(glb2loc_index, &
      &  global_index)
    valid_glb2loc = get_valid_local_index(glb2loc_index, global_index, .TRUE.)
    IF (glb2loc > 0)  THEN
      loc2glb = patch_2D%verts%decomp_info%glb_index(glb2loc)
      blockNo = block_no(glb2loc)
      indexNo = index_no(glb2loc)
      write(0,*) my_proc_id, ": vert: (", patch_2D%verts%vertex(indexNo,blockNo)%lon, &
        & patch_2D%verts%vertex(indexNo, blockNo)%lat, ")", &
        & global_index, glb2loc, valid_glb2loc, loc2glb, owner_local(glb2loc)
      blockNo = block_no(glb2loc)
      indexNo = index_no(glb2loc)
      noOfEdges = patch_2d%verts%num_edges(indexNo, blockNo)
      write(0,*) my_proc_id, ":", global_index, " local edge_of_verts:", &
        & ((patch_2d%verts%edge_blk(indexNo, blockNo, 1:noOfEdges) - 1) * nproma &
        &  + patch_2d%verts%edge_idx(indexNo, blockNo, 1:noOfEdges))
      write(0,*) my_proc_id, ":", global_index, " global edge_of_verts:", &
        & patch_2D%edges%decomp_info%glb_index((patch_2d%verts%edge_blk(indexNo, blockNo, 1:noOfEdges) - 1) * nproma &
        &  + patch_2d%verts%edge_idx(indexNo, blockNo, 1:noOfEdges))
    ENDIF

    global_index = 7290
    glb2loc = get_local_index(glb2loc_index, &
      &  global_index)
    valid_glb2loc = get_valid_local_index(glb2loc_index, global_index, .TRUE.)
    IF (glb2loc > 0)  THEN
      loc2glb = patch_2D%verts%decomp_info%glb_index(glb2loc)
      write(0,*) my_proc_id, ": vert: (", patch_2D%verts%vertex(indexNo,blockNo)%lon, &
        & patch_2D%verts%vertex(indexNo, blockNo)%lat, ")", &
        & global_index, glb2loc, valid_glb2loc, loc2glb, owner_local(glb2loc)
      blockNo = block_no(glb2loc)
      indexNo = index_no(glb2loc)
      noOfEdges = patch_2d%verts%num_edges(indexNo, blockNo)
      write(0,*) my_proc_id, ":", global_index, " local edge_of_verts:", &
        & ((patch_2d%verts%edge_blk(indexNo, blockNo, 1:noOfEdges) - 1) * nproma &
        &  + patch_2d%verts%edge_idx(indexNo, blockNo, 1:noOfEdges))
      write(0,*) my_proc_id, ":", global_index, " global edge_of_verts:", &
        & patch_2D%edges%decomp_info%glb_index((patch_2d%verts%edge_blk(indexNo, blockNo, 1:noOfEdges) - 1) * nproma &
        &  + patch_2d%verts%edge_idx(indexNo, blockNo, 1:noOfEdges))
    ENDIF


    CALL work_mpi_barrier()    

  END SUBROUTINE check_global_indexes
  !----------------------------------------------------

  
END MODULE mo_grid_tools
