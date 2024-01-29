!> @file comin_descrdata_types.F90
!! @brief ComIn descriptive data structures.
!
!  @authors 10/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE comin_descrdata_types
  USE ISO_C_BINDING,         ONLY : C_SIGNED_CHAR
  USE comin_setup_constants, ONLY : wp
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_comin_descrdata_global
  PUBLIC :: t_comin_descrdata_domain
  PUBLIC :: t_comin_descrdata_domain_cells
  PUBLIC :: t_comin_descrdata_domain_edges
  PUBLIC :: t_comin_descrdata_domain_verts
  PUBLIC :: t_comin_descrdata_simulation_interval
  PUBLIC :: comin_glb2loc_index_lookup_fct

#include "comin_global.inc"

  ! ------------------------------------
  ! data types for descriptive data structures
  ! ------------------------------------

  ! Descriptive data from p_patch apart from few exceptions (e.g. lon and lat)
  ! handed to ComIn as pointer, all other descriptive data are copies from
  ! host model, which need to be updated during the model simulation if changed

  !> Global data is invariant wrt the computational grid and never changed
  !! or updated.
  !! @ingroup common
  TYPE :: t_comin_descrdata_global
    ! number of logical domains
    INTEGER               :: n_dom
    ! maximum number of model domains
    INTEGER               :: max_dom
    ! block size
    INTEGER               :: nproma
    ! KIND value (REAL)
    INTEGER               :: wp
    ! block index
    INTEGER               :: min_rlcell_int
    INTEGER               :: min_rlcell
    ! block index
    INTEGER               :: grf_bdywidth_c
    INTEGER               :: grf_bdywidth_e
    ! whether this is a restarted run
    LOGICAL               :: lrestartrun
    ! parameter A of the vertical coordinate (without influence of topography)
    ! index=1,nlev+1
    REAL(wp), ALLOCATABLE :: vct_a(:)
    ! The yac instance id used by the host model
    INTEGER               :: yac_instance_id
  END TYPE t_comin_descrdata_global

  !> Cell information for grid data structures
  TYPE :: t_comin_descrdata_domain_cells
    ! number of local cells
    INTEGER, POINTER :: ncells => NULL()
    ! number of global cells
    INTEGER, POINTER :: ncells_global => NULL()
    ! number of blocks for cells
    INTEGER, POINTER :: nblks => NULL()

    INTEGER, POINTER          :: max_connectivity => NULL()
    ! number of edges connected to cell
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER, CONTIGUOUS          :: num_edges(:,:) => NULL()
    ! lateral boundary distance indices
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER             :: refin_ctrl(:,:) => NULL()
    ! list of start indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell (defined in mo_impl_constants)
    INTEGER, POINTER, CONTIGUOUS          :: start_index(:) => NULL()
    ! list of end indices for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell
    INTEGER, POINTER, CONTIGUOUS          :: end_index(:) => NULL()
    ! list of start block for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell
    INTEGER, POINTER, CONTIGUOUS          :: start_block(:) => NULL()
    ! list of end block for each refin_ctrl level
    ! index1=min_rlcell,max_rlcell
    INTEGER, POINTER, CONTIGUOUS          :: end_block(:) => NULL()
    ! domain ID of child triangles:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER, CONTIGUOUS          :: child_id(:,:) => NULL()
    ! line indices of child triangles:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,4
    INTEGER, POINTER, CONTIGUOUS          :: child_idx(:,:,:) => NULL()
    ! block indices of child triangles:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,4
    INTEGER, POINTER, CONTIGUOUS          :: child_blk(:,:,:) => NULL()
    ! line index of parent triangle:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER, CONTIGUOUS          :: parent_glb_idx(:,:) => NULL()
    ! block index of parent triangle:
    ! index1=1,nproma, index2=1,nblks_c
    INTEGER, POINTER, CONTIGUOUS          :: parent_glb_blk(:,:) => NULL()
    ! line indices and blocks of verts of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, POINTER, CONTIGUOUS          :: vertex_idx(:,:,:) => NULL()
    INTEGER, POINTER, CONTIGUOUS          :: vertex_blk(:,:,:) => NULL()
    ! line indices and blocks of neighbors of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, POINTER, CONTIGUOUS          :: neighbor_blk(:,:,:) => NULL()
    INTEGER, POINTER, CONTIGUOUS          :: neighbor_idx(:,:,:) => NULL()
    ! line indices and blocks of edges of triangle:
    ! index1=1,nproma, index2=1,nblks_c, index3=1,3
    INTEGER, POINTER, CONTIGUOUS          :: edge_idx(:,:,:) => NULL()
    INTEGER, POINTER, CONTIGUOUS          :: edge_blk(:,:,:) => NULL()
    ! longitude & latitude of centers of triangular cells
    ! index1=nproma, index2=1,nblks_c
    REAL(wp), ALLOCATABLE     :: clon(:,:)
    REAL(wp), ALLOCATABLE     :: clat(:,:)
    ! area of triangle
    ! index1=nproma, index2=1,nblks_c
    REAL(wp), POINTER, CONTIGUOUS         :: area(:,:) => NULL()
    ! geometrical height of half levels at cell centre
    ! index1=1,nproma, index2=1,nlev+1, index3=1,nblks_c
    REAL(wp), ALLOCATABLE                 :: hhl(:,:,:)

    ! global cell indices
    INTEGER, POINTER, CONTIGUOUS  :: glb_index(:) => NULL()
    ! Domain decomposition flag:
    ! decomp_domain==0: inner domain, decomp_domain>0: boundary, decomp_domain<0: undefined
    ! For cells:
    ! 0=owned, 1=shared edge with owned, 2=shared vertex with owned
    ! index1=nproma, index2=1,nblks_c
    INTEGER, POINTER, CONTIGUOUS  :: decomp_domain(:,:) => NULL()
  END TYPE t_comin_descrdata_domain_cells

  !> Vertex information for grid data structures
  TYPE :: t_comin_descrdata_domain_verts
    ! number of local verts
    INTEGER, POINTER :: nverts => NULL()
    ! number of global verts
    INTEGER, POINTER :: nverts_global => NULL()
    ! number of blocks for verts
    INTEGER, POINTER :: nblks => NULL()

    ! lateral boundary distance indices
    ! index1=1,nproma, index2=1,nblks_v
    INTEGER, POINTER             :: refin_ctrl(:,:) => NULL()
    ! list of start indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert (defined in mo_impl_constants), index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: start_index(:) => NULL()
    ! list of end indices for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: end_index(:) => NULL()
    ! list of start block for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: start_block(:) => NULL()
    ! list of end block for each refin_ctrl level
    ! index1=min_rlvert,max_rlvert, index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: end_block(:) => NULL()
    ! block indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, POINTER, CONTIGUOUS          :: neighbor_blk(:,:,:) => NULL()
    ! line indices of neighbor vertices:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, POINTER, CONTIGUOUS          :: neighbor_idx(:,:,:) => NULL()
    ! line indices of cells around each vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, POINTER, CONTIGUOUS          :: cell_idx(:,:,:) => NULL()
    ! block indices of cells around each vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, POINTER, CONTIGUOUS          :: cell_blk(:,:,:) => NULL()
    ! line indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, POINTER, CONTIGUOUS          :: edge_idx(:,:,:) => NULL()
    ! block indices of edges around a vertex:
    ! index1=1,nproma, index2=1,nblks_v, index3=1,6
    INTEGER, POINTER, CONTIGUOUS          :: edge_blk(:,:,:) => NULL()
    ! longitude & latitude of vertex:
    ! index1=1,nproma, index2=1,nblks_v
    REAL(wp), ALLOCATABLE     :: vlon(:,:)
    REAL(wp), ALLOCATABLE     :: vlat(:,:)
  END TYPE t_comin_descrdata_domain_verts

  !> Edge information for grid data structures
  TYPE :: t_comin_descrdata_domain_edges
    ! number of local edges
    INTEGER, POINTER :: nedges => NULL()
    ! number of global edges
    INTEGER, POINTER :: nedges_global => NULL()
    ! number of blocks for edges
    INTEGER, POINTER :: nblks => NULL()

    ! lateral boundary distance indices
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, POINTER             :: refin_ctrl(:,:) => NULL()
    ! list of start indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge (defined in mo_impl_constants), index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: start_index(:) => NULL()
    ! list of end indices for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: end_index(:) => NULL()
    ! list of start block for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: start_block(:) => NULL()
    ! list of end block for each refin_ctrl level
    ! index1=min_rledge,max_rledge, index2=n_childdom
    INTEGER, POINTER, CONTIGUOUS          :: end_block(:) => NULL()
    ! domain ID of child edges:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, POINTER, CONTIGUOUS          :: child_id(:,:) => NULL()
    ! line indices of child edges:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, POINTER, CONTIGUOUS          :: child_idx(:,:,:) => NULL()
    ! block indices of child edges:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, POINTER, CONTIGUOUS          :: child_blk(:,:,:) => NULL()
    ! line index of parent edges:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, POINTER, CONTIGUOUS          :: parent_glb_idx(:,:) => NULL()
    ! block index of parent edges:
    ! index1=1,nproma, index2=1,nblks_e
    INTEGER, POINTER, CONTIGUOUS          :: parent_glb_blk(:,:) => NULL()
    ! line indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, POINTER, CONTIGUOUS          :: cell_idx(:,:,:) => NULL()
    ! block indices of adjacent cells:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,2
    INTEGER, POINTER, CONTIGUOUS          :: cell_blk(:,:,:) => NULL()
    ! line indices of edge vertices:
    ! vertex indices 3 and 4 are the non-edge-aligned vertices of cells 1 and 2
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, POINTER, CONTIGUOUS          :: vertex_idx(:,:,:) => NULL()
    ! block indices of edge vertices:
    ! index1=1,nproma, index2=1,nblks_e, index3=1,4
    INTEGER, POINTER, CONTIGUOUS          :: vertex_blk(:,:,:) => NULL()
    ! longitude & latitude of edge midpoint
    ! index=1,nproma, index2=1,nblks_e
    REAL(wp), ALLOCATABLE     :: elon(:,:)
    REAL(wp), ALLOCATABLE     :: elat(:,:)
  END TYPE t_comin_descrdata_domain_edges

  !> Patch grid data structure, gathering information on grids
  !! @ingroup common
  TYPE :: t_comin_descrdata_domain
    ! horizontal grid filename
    CHARACTER(LEN=:), POINTER    :: grid_filename => NULL()
    ! alphanumerical hash of grid
    INTEGER(c_signed_char), POINTER     :: grid_uuid(:) => NULL()
    ! number of grid used (GRIB2 key)
    ! index=1,max_dom
    INTEGER, ALLOCATABLE         :: number_of_grid_used(:)
    ! domain ID of current domain
    INTEGER, POINTER             :: id => NULL()
    ! actual number of child domains
    INTEGER, POINTER             :: n_childdom => NULL()
    ! time at which execution of a (nested) model domain starts
    REAL(wp)                     :: dom_start
    ! time at which execution of a (nested) model domain terminates
    REAL(wp)                     :: dom_end
    ! no. of vertical model levels
    INTEGER, POINTER             :: nlev => NULL()
    ! half level of parent domain that coincides with upper margin of current
    ! domain
    INTEGER, POINTER             :: nshift => NULL()
    ! total shift of model top w.r.t. global domain
    INTEGER, POINTER             :: nshift_total => NULL()

    TYPE(t_comin_descrdata_domain_cells)    :: cells
    TYPE(t_comin_descrdata_domain_verts)    :: verts
    TYPE(t_comin_descrdata_domain_edges)    :: edges
  END TYPE t_comin_descrdata_domain

  !> Simulation status information, sim_current contains current time step
  !! @ingroup common
  TYPE :: t_comin_descrdata_simulation_interval
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: exp_start
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: exp_stop
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: run_start
    CHARACTER(LEN=MAX_DATETIME_STR_LEN) :: run_stop
  END TYPE t_comin_descrdata_simulation_interval

 INTERFACE
    FUNCTION comin_glb2loc_index_lookup_fct(jg, glb) RESULT(loc)
      INTEGER,          INTENT(IN) :: jg
      INTEGER,          INTENT(IN) :: glb
      INTEGER :: loc
    END FUNCTION comin_glb2loc_index_lookup_fct
 END INTERFACE

END MODULE comin_descrdata_types
