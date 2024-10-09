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

! The module contains the definition of the grid and patch data
! structure that is actually used by the model.  It contains abstract
! data types that are copies of those in m_grid, m_hierarchy and
! m_base_geometry.

!! I can't say why xlf takes ages to optimize this file,
!! but since it contains declarations only, it's sufficient to leave
!! it unoptimized
#ifdef __xlC__
@PROCESS NOOPTimize
#endif
MODULE mo_model_domain
  !-------------------------------------------------------------------------
  !
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !
  !-------------------------------------------------------------------------
  !
  !
  !
  USE mo_kind
  USE mo_math_types,              ONLY: t_geographical_coordinates, t_cartesian_coordinates
  USE mo_impl_constants,          ONLY: max_dom, max_phys_dom
  USE mo_communication,           ONLY: t_comm_pattern, t_comm_gather_pattern, &
    &                                   t_scatterPattern, &
    &                                   t_comm_pattern_collection, &
    &                                   t_p_comm_pattern
  USE mo_io_units,                ONLY: filename_max
  USE mo_util_uuid_types,         ONLY: t_uuid
  USE mo_lib_grid_geometry_info,  ONLY: t_grid_geometry_info
  USE mo_decomposition_tools,     ONLY: t_grid_domain_decomp_info
  USE mo_read_netcdf_types,       ONLY: t_distrib_read_data
  USE ppm_distributed_array,      ONLY: dist_mult_array
  USE ppm_extents,                ONLY: extent

  IMPLICIT NONE

  PRIVATE

  ! ! the following abstract data types are taken from mo_grid,

  PUBLIC :: t_patch, t_pre_patch
  PUBLIC :: t_grid_cells, t_pre_grid_cells
  PUBLIC :: t_grid_edges, t_pre_grid_edges
  PUBLIC :: t_grid_vertices
  PUBLIC :: t_phys_patch
  PUBLIC :: t_subset_range, t_subset_range_index, t_subset_indexed

  !PUBLIC :: t_patch_ocean
  PUBLIC :: t_patch_3D
  PUBLIC :: t_patch_vert

  PUBLIC :: t_tangent_vectors

  PUBLIC :: get_startrow_c

  !----------------------------------------------------
  ! TODO Abstraction of the two types of subsets
  !----------------------------------------------------
  !> Defines a subset explicitly using an array of indexes
  TYPE :: t_subset_indexed

    INTEGER, ALLOCATABLE :: block(:)
    INTEGER, ALLOCATABLE :: idx(:)

    INTEGER :: size
    INTEGER :: recommended_stride  ! in case of parallelization/vectorization

    TYPE(t_patch), POINTER :: patch => NULL()
    TYPE(t_patch_3D), POINTER :: patch_3D => NULL()
    INTEGER :: entity_location ! on_cells, on_edges, on_verts

    INTEGER, POINTER :: vertical_levels(:,:) => NULL()  ! if not null, points to the number of vertical levels array

    LOGICAL :: is_in_domain

    CHARACTER(len=32) :: name

  END TYPE t_subset_indexed
  !----------------------------------------------------

  !----------------------------------------------------
  !> Defines a subset in a range (in terms of blocks)
  TYPE :: t_subset_range

    INTEGER :: start_block
    INTEGER :: start_index
    INTEGER :: end_block
    INTEGER :: end_index
    INTEGER :: block_size

    INTEGER :: size

    TYPE(t_patch), POINTER :: patch => NULL()
    INTEGER :: entity_location ! on_cells, on_edges, on_verts

    INTEGER :: max_vertical_levels
    INTEGER, POINTER :: vertical_levels(:,:) => NULL()  ! if not null, points to the number of verticall levels array

    INTEGER :: no_of_holes ! the number of holes in the subset
    LOGICAL :: is_in_domain

    CHARACTER(len=32) :: name

  END TYPE t_subset_range
  !----------------------------------------------------


  !----------------------------------------------------
  !> Defines an index for a subset_range
  TYPE :: t_subset_range_index

    INTEGER, POINTER :: block => NULL() ! the current block in the subset
    INTEGER, POINTER :: index => NULL() ! the current index in the subset

    INTEGER :: start_index ! the current start index within the current block,
    INTEGER :: end_index   ! the current end index within the current block,

    TYPE(t_subset_range), POINTER :: subset_range => NULL()

  END TYPE t_subset_range_index
  !----------------------------------------------------

  ! tangent vector class
  TYPE t_tangent_vectors
    REAL(wp) :: v1
    REAL(wp) :: v2
  END TYPE t_tangent_vectors

  ! -----------------------------------------------------------------------------
  ! grid_cell class - corresponds to triangles
  ! -----------------------------------------------------------------------------
  !
  ! This data type holds the topological and geometrical information on
  ! the grid cells, read from the NetCDF grid file.
  !
  ! * Note on parent-child relation "child_idx/child_blk":
  !
  !          *---------------o--------------*
  !           \__   4     _/  \_    2    __/
  !              \_     _/   3  \_    __/
  !                \__ /__________\ _/
  !                   o            o  
  !                     \__  1  __/
  !                        \_ _/
  !                          *
  !
  !   The ICON model makes an important implicit assumption on the
  !   child ordering:
  !
  !   - cell child #3 is the "interior" child (containing the parent
  !     cell's mass point)
  !   
  TYPE t_grid_cells

    INTEGER :: max_connectivity
    ! number of edges connected to cell
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: num_edges(:,:)

    ! line index of parent triangle:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: parent_loc_idx(:,:)
    INTEGER, ALLOCATABLE :: parent_glb_idx(:,:)
    ! block index of parent triangle:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: parent_loc_blk(:,:)
    INTEGER, ALLOCATABLE :: parent_glb_blk(:,:)
    ! parent child index, number of current cell in parent's child_idx/child_blk:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: pc_idx(:,:)

    ! line indices of child triangles:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,4
    INTEGER, ALLOCATABLE :: child_idx(:,:,:)
    ! block indices of child triangles:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,4
    INTEGER, ALLOCATABLE :: child_blk(:,:,:)
    ! domain ID of child triangles:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: child_id(:,:)
    ! physical domain ID of triangles
    ! (may differ from the "normal" domain ID in case of domain merging):
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: phys_id(:,:)

    ! line indices of triangles next to each cell:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: neighbor_idx(:,:,:)
    ! block indices of triangles next to each cell:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: neighbor_blk(:,:,:)

    ! line indices of edges of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: edge_idx(:,:,:)
    ! block indices of edges of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: edge_blk(:,:,:)

    ! line indices of verts of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: vertex_idx(:,:,:)
    ! block indices of verts of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, ALLOCATABLE :: vertex_blk(:,:,:)

    ! orientation of vectors normal to cell edges:
    ! index1=nproma, index2=1,nblks_c, index3=1,3
    REAL(wp), ALLOCATABLE :: edge_orientation(:,:,:)

    ! cell geometry

    ! longitude & latitude of centers of triangular cells
    ! index1=nproma, index2=1,nblks_c
    TYPE(t_geographical_coordinates), ALLOCATABLE ::  &
      & center(:,:)

    ! area of triangle
    ! index1=nproma, index2=1,nblks_c
    REAL(wp), POINTER :: area(:,:) => NULL()

    ! Coriolis parameter at cell centers
    ! index1=1,nproma, index2=1,nblks_c
    REAL(wp), ALLOCATABLE :: f_c(:,:)

    !----------------------------------
    ! cell geometry auxiliary variables
    ! the cartesian coordinates of the cell centers on the unit sphere
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian_center(:,:) => NULL()
    !----------------------------------

    ! refinement control flag
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell (defined in mo_impl_constants), index2=n_childdom
    INTEGER, ALLOCATABLE :: start_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_index(:) ! revised implementation

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_index(:) ! revised implementation

    ! list of start block for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell, index2=n_childdom
    INTEGER, ALLOCATABLE :: start_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_block(:) ! revised implementation

    ! list of end block for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_block(:) ! revised implementation

    ! information on domain decomposition
    TYPE(t_grid_domain_decomp_info) :: decomp_info

    ! information for distributed read operation
    ! the pointer is needed because mo_read_interface has references this,
    ! without this pointer p_patch would have to be TARGET everywhere the
    ! distributed IO is used
    TYPE(t_distrib_read_data), POINTER :: dist_io_data => NULL()

    ! Please note that the following array is only needed on local parent patches
    ! for storing the corresponding variable from nh_metrics.
    ! It is not allocated/deallocated with the regular patch (de)allocation routines
    ! nor does it need to be dumped/restored or subdivided in patch subdivision.
    REAL(wp), ALLOCATABLE :: ddqz_z_full(:,:,:)

    ! define basic subsets
    TYPE(t_subset_range) :: ALL          ! these are the all valid entities, including all valid halos
    TYPE(t_subset_range) :: owned         ! these are the owned entities
    TYPE(t_subset_range) :: owned_no_boundary ! these are the owned entities without boundary interpolation zone
    TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities, for cells = owned
    TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
    TYPE(t_subset_range) :: not_in_domain ! for cells = not_owned
    TYPE(t_subset_range) :: one_edge_in_domain ! cells with exactly one edge in domain. these are always halo cells

    INTEGER :: dummy_cell_block, dummy_cell_index     ! in case of existing a dummy cell for closing the
      ! local domain, these point to the dummy cell. This cell does not actually "exist"!
      ! =0 if no dummy cell

  END TYPE t_grid_cells


  ! index of distributed sub-arrays in t_pre_grid_cells%dist
  INTEGER, PUBLIC, PARAMETER :: &
       c_num_edges = 8, &    ! number of edges connected to cell
       c_parent = 2, &    ! index of parent triangle:
       ! parent child index, number of current cell in parent's child_idx/child_blk:
       ! indices of child triangles:
       ! index2=1,4
       c_child = 3, &
       ! physical domain ID of triangles
       ! (may differ from the "normal" domain ID in case of domain merging):
       c_phys_id = 4, &
       ! indices of triangles next to each cell:
       ! index2=1,3
       c_neighbor = 5, &
       ! indices of edges of triangle:
       ! index2=1,3
       c_edge = 6, &
       ! indices of verts of triangle:
       ! index2=1,3
       c_vertex = 7, &
       ! cell geometry
       ! longitude & latitude of centers of triangular cells
       c_center = 1, &
       ! refinement control flag
       c_refin_ctrl = 9

  TYPE t_pre_grid_cells

    ! extents of the local chunk of the distributed arrays
    TYPE(extent) :: local_chunk(1,1)

    INTEGER :: max_connectivity
    TYPE(dist_mult_array) :: dist

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell (defined in mo_impl_constants)
    INTEGER, ALLOCATABLE :: start(:)

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell
    INTEGER, ALLOCATABLE :: end(:)

  END TYPE t_pre_grid_cells

  ! !grid_edge class

  TYPE t_grid_edges

    ! line index of parent edge:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: parent_loc_idx(:,:)
    INTEGER, ALLOCATABLE :: parent_glb_idx(:,:)
    ! block index of parent edge:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: parent_loc_blk(:,:)
    INTEGER, ALLOCATABLE :: parent_glb_blk(:,:)
    ! parent child index, number of current edge in parent's child_idx/child_blk:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: pc_idx(:,:)

    ! line indices of child edges:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: child_idx(:,:,:)
    ! block indices of child edges:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: child_blk(:,:,:)
    ! domain ID of child edges:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: child_id(:,:)
    ! physical domain ID of edges
    ! (may differ from the "normal" domain ID in case of domain merging):
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: phys_id(:,:)

    ! line indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, ALLOCATABLE :: cell_idx(:,:,:)
    ! block indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, ALLOCATABLE :: cell_blk(:,:,:)

    ! line indices of edge vertices:
    ! vertex indices 3 and 4 are the non-edge-aligned vertices of cells 1 and 2
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: vertex_idx(:,:,:)
    ! block indices of edge vertices:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: vertex_blk(:,:,:)

    ! =1 if vector product of vector from vertex1 to vertex 2 (v2-v1) by vector
    ! from cell c1 to cell c2 (c2-c1) goes outside the sphere
    ! =-1 if vector product ...       goes inside  the sphere
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: tangent_orientation(:,:)

    ! line indices of the  of the quadrilateral formed by two adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: quad_idx(:,:,:)
    ! block indices of the  of the quadrilateral formed by two adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, ALLOCATABLE :: quad_blk(:,:,:)

    ! orientation of vectors normal to quad cell edges:
    ! index1=nproma, index2=1,nblks_e, index3=1,4
    REAL(wp), ALLOCATABLE :: quad_orientation(:,:,:)

    ! line indices of cells which share the two vertices
    ! of a given edge. These cells are subdivided into two classes:
    ! Cells which are neighbors of edge neighbor 1 and cells which are
    ! neighbors of edge neighbor 2. These cells are then numbered
    ! according to the number of the edge-vertex they share.
    ! index1=nproma, index2=1,nblks_e, index3=1,2 (cell), index4=1,2 (vert)
    INTEGER, ALLOCATABLE :: butterfly_idx(:,:,:,:)
    ! block indices of cells which share the two vertices
    ! of a given edge. These cells are subdivided into two classes:
    ! Cells which are neighbors of edge neighbor 1 and cells which are
    ! neighbors of edge neighbor 2. These cells are then numbered
    ! according to the number of the edge-vertex they share.
    ! index1=nproma, index2=1,nblks_e, index3=1,2 (cell), index4=1,2 (vert)
    INTEGER, ALLOCATABLE :: butterfly_blk(:,:,:,:)


    !-------------------------------------------------
    ! edges geometry

    ! longitude & latitude of edge midpoint
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_geographical_coordinates), ALLOCATABLE ::  &
      & center(:,:)

    ! normal to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & primal_normal(:,:)

    ! Cartesian normal to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & primal_cart_normal(:,:) => NULL()

    ! Cartesian dual to triangle edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & dual_cart_normal(:,:) => NULL()

    ! normal to hexagon/pentagon edge
    ! index=1,nproma, index2=1,nblks_e
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & dual_normal(:,:)

    ! normal to triangle edge, projected to the location of the neighbor cells
    ! index=1,nproma, index2=1,nblks_e, index3=1,2
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & primal_normal_cell(:,:,:)

    ! tangent to triangle edge, projected to the location of the neighbor cells
    ! index=1,nproma, index2=1,nblks_e, index3=1,2
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & dual_normal_cell(:,:,:)

    ! normal to triangle edge, projected to the location of the vertices
    ! index=1,nproma, index2=1,nblks_e, index3=1,4
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & primal_normal_vert(:,:,:)

    ! tangent to triangle edge, projected to the location of the vertices
    ! index=1,nproma, index2=1,nblks_e, index3=1,4
    TYPE(t_tangent_vectors), ALLOCATABLE ::  &
      & dual_normal_vert(:,:,:)

    ! length of triangle edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), POINTER :: primal_edge_length(:,:) => NULL()

    ! inverse length of triangle edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: inv_primal_edge_length(:,:)

    ! length of hexagon/pentagon edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: dual_edge_length(:,:)

    ! inverse length of hexagon/pentagon edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: inv_dual_edge_length(:,:)

    ! length of edge midpoint to vertex
    ! index=1,nproma, index2=1,nblks_e,2
    REAL(wp), ALLOCATABLE :: edge_vert_length(:,:,:)

    ! length of edge midpoint to cell center
    ! index=1,nproma, index2=1,nblks_e,2
    REAL(wp), ALLOCATABLE :: edge_cell_length(:,:,:)

    ! inverse distance between outer vertices of adjacent cells
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: inv_vert_vert_length(:,:)

    ! area of the quadrilateral given by the adjacent vertices and centers
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: area_edge(:,:)

    ! area of the quadrilateral formed by two cells adjacent to the edge
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: quad_area(:,:)

    ! edge geometry auxiliary variables
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian_center(:,:) => NULL()

    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian_dual_middle(:,:) => NULL()

    ! Coriolis parameter at cell edges
    ! index1=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE :: f_e(:,:)

    ! components of horizontal Coriolis parameter (vector)
    ! index1=1,nproma, index2=1,nblks_e
    ! NOTE: The horizontal Coriolis parameter components 
    ! are only available on edges!
    ! normal component 
    REAL(wp), ALLOCATABLE :: fn_e(:,:)
    ! tangential component
    REAL(wp), ALLOCATABLE :: ft_e(:,:)

    ! refinement control flag
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

    ! list of start indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge (defined in mo_impl_constants), index2=n_childdom
    INTEGER, ALLOCATABLE :: start_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_index(:) ! revised implementation

    ! list of end indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_index(:) ! revised implementation

    ! list of start block for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, ALLOCATABLE :: start_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_block(:) ! revised implementation

    ! list of end block for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_block(:) ! revised implementation

    ! information on domain decomposition
    TYPE(t_grid_domain_decomp_info) :: decomp_info

    ! information for distributed read operation
    ! the pointer is needed because mo_read_interface has references this,
    ! without this pointer p_patch would have to be TARGET everywhere the
    ! distributed IO is used
    TYPE(t_distrib_read_data), POINTER :: dist_io_data => NULL()

    ! define basic subsets
    TYPE(t_subset_range) :: ALL          ! these are the all valid entities, including all valid halos
    TYPE(t_subset_range) :: owned         ! these are the owned entities
    TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities,
    TYPE(t_subset_range) :: gradIsCalculable     ! these are edges for which the grad can be calcluated
    ! this includes all edges of the domain cells
    TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
    TYPE(t_subset_range) :: not_in_domain ! all - in_domain

  END TYPE t_grid_edges

  ! index of distributed sub-arrays in t_pre_grid_edges%dist
  INTEGER, PUBLIC, PARAMETER :: &
       ! index of parent edge:
       e_parent = 1, &
       ! indices of child edges:
       ! index2=1,4
       e_child = 2, &
       ! indices of adjacent cells:
       ! index2=1,2
       e_cell = 3, &
       ! refinement control flag
       e_refin_ctrl = 4

  TYPE t_pre_grid_edges

    ! extents of the local chunk of the distributed arrays
    TYPE(extent) :: local_chunk(1,1)

    TYPE(dist_mult_array) :: dist

    ! list of start indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge (defined in mo_impl_constants)
    INTEGER, ALLOCATABLE :: start(:)

    ! list of end indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge
    INTEGER, ALLOCATABLE :: end(:)

  END TYPE t_pre_grid_edges

  !  !grid_vertices class

  TYPE t_grid_vertices

    INTEGER :: max_connectivity
    ! physical domain ID of verts
    ! (may differ from the "normal" domain ID in case of domain merging):
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, ALLOCATABLE :: phys_id(:,:)

    ! line indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: neighbor_idx(:,:,:)
    ! block indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: neighbor_blk(:,:,:)

    ! line indices of cells around each vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: cell_idx(:,:,:)
    ! block indices of cells around each vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: cell_blk(:,:,:)

    ! line indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: edge_idx(:,:,:)
    ! block indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, ALLOCATABLE :: edge_blk(:,:,:)

    ! (xx)
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    REAL(wp), ALLOCATABLE :: edge_orientation(:,:,:)

    ! number of edges connected to vertex
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, ALLOCATABLE :: num_edges(:,:)

    ! longitude & latitude of vertex:
    ! index1=1,nproma, index2=1,nblks_v
    TYPE(t_geographical_coordinates), ALLOCATABLE :: vertex(:,:)

    ! area of hexagon/pentagon of which vertex is center:
    ! index1=1,nproma, index2=1,nblks_v
    REAL(wp), POINTER :: dual_area(:,:) => NULL()

    ! Coriolis parameter at cell vertices
    ! index1=1,nproma, index2=1,nblks_v
    REAL(wp), ALLOCATABLE :: f_v(:,:)

    ! vertex geometry auxiliary variables
    TYPE(t_cartesian_coordinates), POINTER ::  &
      & cartesian(:,:) => NULL()

    ! refinement control flag
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, ALLOCATABLE :: refin_ctrl(:,:)

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert (defined in mo_impl_constants), index2=n_childdom
    INTEGER, ALLOCATABLE :: start_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_index(:) ! revised implementation

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_idx(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_index(:) ! revised implementation

    ! list of start block for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, ALLOCATABLE :: start_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: start_block(:) ! revised implementation

    ! list of end block for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, ALLOCATABLE :: end_blk(:,:) ! to be removed soon
    INTEGER, ALLOCATABLE :: end_block(:) ! revised implementation

    ! information on domain decomposition
    TYPE(t_grid_domain_decomp_info) :: decomp_info

    ! information for distributed read operation
    ! the pointer is needed because mo_read_interface has references this,
    ! without this pointer p_patch would have to be TARGET everywhere the
    ! distributed IO is used
    TYPE(t_distrib_read_data), POINTER :: dist_io_data => NULL()

    ! define basic subsets
    TYPE(t_subset_range) :: ALL          ! these are the all valid entities, including all valid halos
    TYPE(t_subset_range) :: owned         ! these are the owned entities
    TYPE(t_subset_range) :: in_domain     ! these are the computation domain entities,
    ! this includes all edges of the domain cells
    TYPE(t_subset_range) :: not_owned     ! these are all the halo entities
    TYPE(t_subset_range) :: not_in_domain ! all - in_domain

  END TYPE t_grid_vertices

  ! index of distributed sub-arrays in t_pre_grid_vertices%dist
  INTEGER, PUBLIC, PARAMETER :: &
       ! line indices of cells around each vertex:
       ! index2=1,6
       v_cell = 3, &
       ! number of edges connected to vertex
       v_num_edges = 2, &
       ! longitude & latitude of vertex:
       ! index2=1,2
       v_vertex = 1, &
       ! refinement control flag
       v_refin_ctrl = 4

  TYPE t_pre_grid_vertices

    ! extents of the local chunk of the distributed arrays
    TYPE(extent) :: local_chunk(1,1)

    INTEGER :: max_connectivity

    TYPE(dist_mult_array) :: dist

    ! list of start indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert (defined in mo_impl_constants)
    INTEGER, ALLOCATABLE :: start(:)

    ! list of end indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert
    INTEGER, ALLOCATABLE :: end(:)

  END TYPE t_pre_grid_vertices

  ! !patch class

  TYPE t_patch

    !
    ! !  level in grid hierarchy on which patch lives
    !
    CHARACTER(LEN=filename_max) :: grid_filename, grid_filename_grfinfo
    !
    ! uuid of grid
    TYPE(t_uuid) :: grid_uuid
    !
    ! grid level
    INTEGER :: level
    !
    ! nest level = grid level - start level
    INTEGER :: nest_level
    !
    ! domain ID of current domain
    INTEGER :: id
    !
    ! indicator if current model domain is active
    LOGICAL :: ldom_active
    !
    ! indicator if current model domain is owned by a PE (relevant in case of processor splitting)
    LOGICAL :: domain_is_owned

    !-------------------------------------
    !> The grid domain geometry parameters
    ! cell type =3 or 6
    ! INTEGER :: cell_type     ! included in geometry_info

!     INTEGER :: geometry_type ! included in geometry_info

    TYPE(t_grid_geometry_info) :: geometry_info
    !-------------------------------------
    INTEGER :: boundary_depth_index  ! when is limited area grid, this is the  number of boundary levels based
                                     ! on the edge-connected cells (the fisrt level are cells that have
                                     ! at least one boundary edge, the next level are cells shering an edge with
                                     ! level 1 cells, etc

    !
    ! domain ID of parent domain
    INTEGER :: parent_id
    !
    ! child domain index of current domain as seen from parent domain
    ! In other words: I am the nth child of my parents (n=parent_child_index)
    INTEGER :: parent_child_index
    !
    ! list of child domain ID's
    INTEGER :: child_id(max_dom)
    !
    ! actual number of child domains
    INTEGER :: n_childdom
    !
    ! total number of child domains in the calling tree (over all nest levels)
    INTEGER :: n_chd_total
    !
    ! corresponding list of child domain ID's
    INTEGER :: child_id_list(max_dom)
    !
    ! maximum number of child domains
    INTEGER :: max_childdom
    !

    ! total number of locally allocated cells, edges and vertices
    INTEGER :: n_patch_cells
    INTEGER :: n_patch_edges
    INTEGER :: n_patch_verts

    ! total number of exist cells
!    INTEGER :: n_exist_cells
!    INTEGER :: n_cell_blocks
!    INTEGER :: last_cell_block_size

    ! in case of a dummy cell, we keep the blk and index of it
!    INTEGER :: dummy_cell_blk
!    INTEGER :: dummy_cell_idx


    !
    ! ! number of cells, edges and vertices in the global patch
    INTEGER :: n_patch_cells_g
    INTEGER :: n_patch_edges_g
    INTEGER :: n_patch_verts_g
    !
    ! ! values for the blocking
    !
    INTEGER :: alloc_cell_blocks  ! number of allocated cell blocks
    ! number of blocks and chunk length in last block
    ! ... for the cells
    INTEGER :: nblks_c
    INTEGER :: npromz_c
    ! ... for the edges
    INTEGER :: nblks_e
    INTEGER :: npromz_e
    ! ... for the vertices
    INTEGER :: nblks_v
    INTEGER :: npromz_v
    !
    ! ! vertical full and half levels
    !
    ! number of full and half levels
    INTEGER :: nlev
    INTEGER :: nlevp1
    !
    ! half level of parent domain (jg-1) that coincides
    ! with the upper margin of the current domain jg
    INTEGER :: nshift
    !
    ! total shift of model top with respect to global domain
    INTEGER :: nshift_total
    !
    ! the same information seen from the parent level (duplication needed to simplify flow control)
    INTEGER :: nshift_child


    !
    ! ! Mask and bathymetry for ocean patch
    !  TYPE(t_patch_ocean) ::  &
    !    &  patch_oce

    !
    ! ! grid information on the patch
    !
    TYPE(t_grid_cells) ::  &
      & cells
    TYPE(t_grid_edges) ::  &
      & edges
    TYPE(t_grid_vertices) ::  &
      & verts

    !
    ! communication patterns for parallelization
    !
    ! Boundary exchange within patches, defined on regular patches and local parents
    CLASS(t_comm_pattern), POINTER :: comm_pat_c
    CLASS(t_comm_pattern), POINTER :: comm_pat_c1 ! reduced communication pattern, only level-1 halo cells
    CLASS(t_comm_pattern), POINTER :: comm_pat_e
    CLASS(t_comm_pattern), POINTER :: comm_pat_v

    ! Interpolation for grid refinement, defined only on regular patches
    CLASS(t_comm_pattern), POINTER :: comm_pat_interpolation_c
    CLASS(t_comm_pattern_collection), POINTER :: comm_pat_coll_interpol_vec_grf
    CLASS(t_comm_pattern_collection), POINTER :: comm_pat_coll_interpol_scal_grf
    CLASS(t_comm_pattern_collection), POINTER :: comm_pat_coll_interpol_vec_ubc
    CLASS(t_comm_pattern_collection), POINTER :: comm_pat_coll_interpol_scal_ubc

    ! Gather complete patch to proc 0
    ! Useful only for regular patches (defined but unused on local parents)
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_c
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_e
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_v

    ! Scatter complete patch from proc 0
    ! Useful only for regular patches (defined but unused on local parents)
    CLASS(t_scatterPattern), POINTER :: comm_pat_scatter_c => NULL()
    CLASS(t_scatterPattern), POINTER :: comm_pat_scatter_e => NULL()
    CLASS(t_scatterPattern), POINTER :: comm_pat_scatter_v => NULL()

    !> if p_test_run, points to array containing
    !! comm_pat_work2test(1) pattern to send cell data to test processes
    !! comm_pat_work2test(2) pattern to send vertex data to test processes
    !! comm_pat_work2test(3) pattern to send edge data to test processes
    TYPE(t_p_comm_pattern) :: comm_pat_work2test(3)

    ! Communication between local parent and its global counterpart,
    ! defined only on local parents.
    ! Please note that these communicate between
    ! p_patch_local_parent(jg) and p_patch(p_patch(jg)%parent_id)
    CLASS(t_comm_pattern), POINTER :: comm_pat_glb_to_loc_c
    CLASS(t_comm_pattern), POINTER :: comm_pat_glb_to_loc_e
    CLASS(t_comm_pattern), POINTER :: comm_pat_loc_to_glb_c_fbk
    CLASS(t_comm_pattern), POINTER :: comm_pat_loc_to_glb_e_fbk

    ! Halo comm patterns for the icon_commm_lib
    INTEGER ::  sync_cells_not_in_domain
    INTEGER ::  sync_cells_not_owned       ! = cells_not_in_domain
    INTEGER ::  sync_cells_one_edge_in_domain
    INTEGER ::  sync_edges_not_owned
    INTEGER ::  sync_edges_not_in_domain
    INTEGER ::  sync_verts_not_owned
    INTEGER ::  sync_verts_not_in_domain

    ! basic parallelization info for this patch
    LOGICAL :: compute_is_parallel  ! if true more than 1 processes are used
    LOGICAL :: is_in_parallel_test  ! if true there's one process that runs seq and compares the results to the parallel runs

    ! basic communicators concerning this patch
    !> the work universe of this patch, halo exchange, global sum,
    ! etc, takes place in this universe


    INTEGER :: work_communicator
    !> the work universe plus the sequential test process,
    ! in case it exists (otherwise the same as work_communicator)
    INTEGER :: parallel_test_communicator

    ! MPI communicator for this patch
    INTEGER :: comm
    !
    ! my index among the processors on this patch
    INTEGER :: rank
    !
    ! total number of processors on this patch
    INTEGER :: n_proc
    !
    ! global id of processor with rank 0 (within the working set p_comm_work)
    INTEGER :: proc0

  END TYPE t_patch
  !-----------------------------------------------------------------------------------
  ! !prepare patch class

  TYPE t_pre_patch

    !
    ! !  level in grid hierarchy on which patch lives
    !
    CHARACTER(LEN=filename_max) :: grid_filename, grid_filename_grfinfo
    !
    ! uuid of grid
    TYPE(t_uuid) :: grid_uuid
    !
    ! grid level
    INTEGER :: level
    !
    ! nest level = grid level - start level
    INTEGER :: nest_level
    !
    ! domain ID of current domain
    INTEGER :: id

    !-------------------------------------
    !> The grid domain geometry parameters
    ! cell type =3 or 6
    INTEGER :: cell_type

    !
    ! domain ID of parent domain
    INTEGER :: parent_id
    !
    ! child domain index of current domain as seen from parent domain
    ! In other words: I am the nth child of my parents (n=parent_child_index)
    INTEGER :: parent_child_index
    !
    ! list of child domain ID's
    INTEGER :: child_id(max_dom)
    !
    ! actual number of child domains
    INTEGER :: n_childdom
    !
    ! total number of child domains in the calling tree (over all nest levels)
    INTEGER :: n_chd_total
    !
    ! corresponding list of child domain ID's
    INTEGER :: child_id_list(max_dom)
    !
    ! maximum number of child domains
    INTEGER :: max_childdom
    !
    ! ! number of cells, edges and vertices in the global patch
    INTEGER :: n_patch_cells_g
    INTEGER :: n_patch_edges_g
    INTEGER :: n_patch_verts_g
    !
    ! ! values for the blocking
    !
    INTEGER :: alloc_cell_blocks  ! number of allocated cell blocks
    ! number of blocks and chunk length in last block
    ! ... for the cells
    INTEGER :: nblks_c
    INTEGER :: npromz_c
    ! ... for the edges
    INTEGER :: nblks_e
    INTEGER :: npromz_e
    ! ... for the vertices
    INTEGER :: nblks_v
    INTEGER :: npromz_v
    !
    ! ! vertical full and half levels
    !
    ! number of full and half levels
    INTEGER :: nlev
    INTEGER :: nlevp1
    !
    ! half level of parent domain (jg-1) that coincides
    ! with the upper margin of the current domain jg
    INTEGER :: nshift
    !
    ! total shift of model top with respect to global domain
    INTEGER :: nshift_total
    !
    ! the same information seen from the parent level (duplication needed to simplify flow control)
    INTEGER :: nshift_child
    !
    ! ! grid information on the patch
    !
    TYPE(t_pre_grid_cells) :: cells
    TYPE(t_pre_grid_edges) :: edges
    TYPE(t_pre_grid_vertices) :: verts
    !
    ! Most of the data of the cells, edges and verts in t_pre_patch is stored in
    ! dist_mult_array. For performance reasons the process in p_comm_work are
    ! split into one or multiple contiguous chunks. Together the processes of
    ! each chunk has a complete copy of the cells, edges and verts data. The
    ! number of chunks is set by the name list parameter num_dist_array_replicas
    !
    ! ! communicator containing all process of the process chunk
    !
    INTEGER :: dist_array_comm
    !
    ! ! process ranks (+1) of the chunk the local process is a part of
    !
    TYPE(extent) :: dist_array_pes

  END TYPE t_pre_patch

  !-----------------------------------------------------------------------------------
  ! Description of physical patches
  TYPE t_phys_patch
    !
    ! domain ID of associated logical domain
    INTEGER :: logical_id
    !
    ! total number of cells, edges and vertices (always global)
    INTEGER :: n_patch_cells
    INTEGER :: n_patch_edges
    INTEGER :: n_patch_verts

    ! Gather the physical patch to proc 0
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_c
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_e
    TYPE(t_comm_gather_pattern) :: comm_pat_gather_v

    ! Scatter the physical patch from proc 0
    ! Currently there is no need for this, so we avoid the overhead of creating it.
    ! Add creation to mo_complete_subdivision::setup_phys_patches and mo_complete_subdivision::setup_phys_patches_cve if you need
    ! this.
!   CLASS(t_scatterPattern), POINTER :: comm_pat_scatter_c => NULL()
!   CLASS(t_scatterPattern), POINTER :: comm_pat_scatter_e => NULL()
!   CLASS(t_scatterPattern), POINTER :: comm_pat_scatter_v => NULL()

  END TYPE t_phys_patch

  TYPE(t_patch), PUBLIC, TARGET, ALLOCATABLE :: p_patch(:)


  ! Definition of local parent patches
  ! For any given patch p_patch(jg) and jgp = p_patch(jg)%parent_id,
  ! p_patch_local_parent(jg) has the same resolution as p_patch(jgp)
  ! but it covers only the area of p_patch(jgp) which is covered by its child p_patch(jg)
  ! and it is divided in the same manner as p_patch(jg).
  ! Please note that p_patch_local_parent(1) is undefined if n_dom_start = 1

  TYPE(t_patch), PUBLIC, TARGET, ALLOCATABLE :: p_patch_local_parent(:)


  ! Please note: There is currently no means of determining the number
  ! of physical patches until they are actually assembled
  ! since this number is missing in the input file.
  ! Therefore p_phys_patch is declared with its maximum size.
  ! This shouldn't hurt since since TYPE t_phys_patch is small.

  TYPE(t_phys_patch), PUBLIC, TARGET :: p_phys_patch(max_phys_dom)

  !--------------------------------------------------------------------

  TYPE t_patch_vert

    !! The ocean uses z-coordinates in meters in the vertical.
    !! The following data are required:
    !!
    !! n_zlev: number of z-coordinate surfaces
    !! n_zlvp: number of intermediate levels (+1)
    !! n_zlvm: number of z-coordinate distances (-1)
    INTEGER :: n_zlev, n_zlvp, n_zlvm

    !! del_zlev_m: thickness (height) of elemental prism, defined as the
    !!             distance between top and bottom of elemental prism,
    !!             i.e. the distance between two intermediate z-coordinate
    !!             surfaces. These data are provided by the user, all other
    !!             vertical information is calculated from this array of
    !!             thicknesses.  Dimension: n_zlev
    REAL(wp), ALLOCATABLE :: del_zlev_m(:)

    !! the inverse of the above
    REAL(wp), ALLOCATABLE :: inv_del_zlev_m(:)

    !! zlev_m    : position of the vertical cell centers, i.e. below zero surface;
    !!             Numbering starts from surface and increases downwards to bottom.
    !!             Dimension: n_zlev
    !!             At these surfaces the horizontal velocities, vorticity, divergence
    !!             and scalar variables are evaluated.
    REAL(wp), ALLOCATABLE :: zlev_m(:)

    !! zlev_i    : vertical position of the UPPER BORDER of the vertical cell
    !!             i.e. the position of the top of elemental prisms.
    !!             Position of first surface is 0.
    !!             Dimension: n_zlvp = n_zlev + 1
    !!             The vertical velocities are evaluated at such surfaces.
    REAL(wp), ALLOCATABLE :: zlev_i(:)

    !! del_zlev_i: distance between two z-coordinate surfaces. The first is
    !!             the distance from the ocean surface = zlev_m(1)
    !!             Dimension: n_zlev
    REAL(wp), ALLOCATABLE :: del_zlev_i(:)

    ! To simplify the acess to the required information within these loops
    ! we store an cell and edge based version of the deepest ocean layer
    ! in column. dolic_e(edge1) and dolic_c(cell1) are identical if 'edge1'
    ! is one of the edges of 'cell1'.
    ! If the ocean bottom is flat dolic_c and dolic_e are identical and equal
    ! to the number of z-coodinate surfaces.
    !
    INTEGER, POINTER :: dolic_c(:,:) => NULL()    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER :: dolic_e(:,:)  => NULL()   ! index1=1,nproma, index2=1,nblks_e
    INTEGER, POINTER :: vertex_bottomLevel(:,:) => NULL()

    REAL(wp), POINTER :: ocean_area   (:) => NULL()  ! global ocean area for each vertical level
    REAL(wp), POINTER :: ocean_volume (:) => NULL()  ! global ocean volume for each vertical level

    REAL(wp), POINTER ::                &
      &  prism_thick_c(:,:,:),          & ! individual prism thickness at cells. Unit [m].
                                          ! This array includes the free surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_volume(:,:,:),           & ! individual prism volume (cells). Unit [m^3].
                                          ! This array includes the free surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_thick_e(:,:,:),          & ! individual prism thickness at edges. Unit [m].
                                          ! This array includes the free surface. dimension: (nproma, n_zlev, nblks_e)
      &  prism_thick_flat_sfc_c(:,:,:) ,& ! individual fluid column thickness at cells. Unit [m].
                                          ! This array assumes a flat surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_thick_flat_sfc_e(:,:,:) ,& ! individual fluid column thickness at edges. Unit [m].
                                          ! This array assumes a flat surface. dimension: (nproma, n_zlev, nblks_c)
      &  prism_center_dist_c(:,:,:),    & ! distance between prism centers at cells. Unit [m].
                                          ! dimension: (nproma,n_zlev, nblks_c)
      &  inv_prism_thick_c(:,:,:) ,     & ! inverse individual fluid column thickness at cells. Unit [m].
                                          ! This array assumes a flat surface dimension: (nproma, n_zlev, nblks_c)
      &  inv_prism_thick_e(:,:,:) ,     & ! inverse individual fluid column thickness at edges. Unit [m].
                                          ! This array assumes a flat surface. dimension: (nproma, n_zlev, nblks_e)
      &  inv_prism_center_dist_c(:,:,:),& ! inverse vertical distance between prism centers at cells. Unit [m].
                                          ! This array assumes a flat surface  dimension: (nproma, n_zlev, nblks_c)
      &  inv_prism_center_dist_e(:,:,:),& ! inverse vertical distance between prism centers at edges. Unit [m].
                                          ! This array assumes a flat surface dimension: (nproma, n_zlev, nblks_e)
      &  depth_CellMiddle(:,:,:),      & ! depth of the middle of the prism, update according to the current h
      &  depth_CellInterface(:,:,:)    & ! depths at the interface (size is levels + 1)
      &  => NULL()
    REAL(wp), POINTER :: invConstantPrismThickness(:,:,:) => NULL()
    !! the vertical distance between the prism centers, dim = n_zlev+1
    !! constantPrismCenters_zDistance(2) = vertical distance between prisms at 1 and 2 levels
    REAL(wp), POINTER :: constantPrismCenters_Zdistance(:,:,:) => NULL()
    !!the inverse of the above
    REAL(wp), POINTER :: constantPrismCenters_invZdistance(:,:,:) => NULL()
    
  END TYPE t_patch_vert



  TYPE t_patch_3D

    TYPE(t_patch),     POINTER :: p_patch_2D(:) => NULL()
    TYPE(t_patch_vert),POINTER :: p_patch_1D(:) => NULL()

    ! land-sea-mask for ocean has 3 dimensions (the 2nd is the number of
    ! vertical levels)
    ! sea=-2, sea_boundary=-1, boundary (edges only)=0, land_boundary=1, land=2
    !
    ! land-sea-mask for cell centers
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_c
    INTEGER, POINTER :: lsm_c(:,:,:) => NULL()

    ! land-sea-mask for edges
    ! index1=1,nproma, index2=1,n_zlev, index3=1,nblks_e
    INTEGER, POINTER :: lsm_e(:,:,:) => NULL()

    ! surface land-sea-mask for cells, edges, vertices
    INTEGER, POINTER :: surface_cell_sea_land_mask(:,:) => NULL()
    INTEGER, POINTER :: surface_edge_sea_land_mask(:,:) => NULL()
    INTEGER, POINTER :: surface_vertex_sea_land_mask(:,:) => NULL()

    ! To simply set land points to zero we store additional 3-dim wet points
    ! dimensions as in lsm_oce:
    REAL(wp), POINTER :: wet_c(:,:,:) => NULL()  ! cells
    REAL(wp), POINTER :: wet_e(:,:,:) => NULL()  ! edges
    ! For calculation of global sum and area including lsm the halo must be set to zero:
    REAL(wp), POINTER :: wet_halo_zero_c(:,:,:) => NULL()  !  cells
    REAL(wp), POINTER :: wet_halo_zero_e(:,:,:) => NULL()  !  edges

    ! For diagnosis like stream functions and area calculations we add surface arrays
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER,  POINTER :: basin_c(:,:) => NULL()  ! basin information Atlantic/Indian/Pacific
    INTEGER,  POINTER :: regio_c(:,:) => NULL()  ! area information like tropical Atlantic etc.
    REAL(wp), POINTER :: bottom_thick_c(:,:) => NULL()  ! individual bottom prism thickness at cells. Unit [m].
    REAL(wp), POINTER :: bottom_thick_e(:,:) => NULL()  ! individual bottom prism thickness at edges. Unit [m].
    REAL(wp), POINTER :: column_thick_c(:,:) => NULL()  ! individual column thickness at cells, no elevation. Unit [m].
    REAL(wp), POINTER :: column_thick_e(:,:) => NULL()  ! individual column thickness at edges, no elevation. Unit [m].

  END TYPE t_patch_3D
  !--------------------------------------------------------------------


CONTAINS

  ! Determine the start row for cells from a given edge start row.
  !
  ! The cell start row defines the row at which cell-based computations start.
  !
  ! For a given edge start row, the cell start row is computed as the minimum of the two
  ! cell rows adjacent to the given edge row. Hence, the cell row which is closer
  ! to the boundary is chosen.
  !
  !
  ELEMENTAL FUNCTION get_startrow_c (startrow_e)  RESULT (startrow_c)

    INTEGER, INTENT(IN) :: startrow_e  !< row at which edge-based computations start
    INTEGER :: startrow_c

    startrow_c = MAX(1,(startrow_e - MOD(startrow_e,2)) / 2)

  END FUNCTION get_startrow_c

END MODULE mo_model_domain

