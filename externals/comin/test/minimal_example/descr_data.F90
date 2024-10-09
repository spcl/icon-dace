!> Data structures and procedures resembling ICON descriptive data structures.
!
!  @authors 08/2021 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
!
MODULE descr_data

  USE, INTRINSIC :: iso_c_binding, ONLY: c_signed_char
  USE comin_host_interface,   ONLY : t_comin_descrdata_global, comin_descrdata_set_global,     &
    &                                t_comin_descrdata_domain, comin_descrdata_set_domain,     &
    &                                t_comin_descrdata_simulation_interval, comin_descrdata_set_simulation_interval, &
    &                                comin_current_set_datetime, comin_descrdata_set_fct_glb2loc_cell, &
    &                                comin_descrdata_set_timesteplength
  USE mo_utilities,           ONLY : wp,idx_no, blk_no, finish, t_geographical_coordinates, &
    &                                read_netcdf_scalar, filename_max, MAX_DATETIME_STR_LEN
  USE mpi,                    ONLY : MPI_Allreduce, MPI_Comm_rank,&
    &                                MPI_Comm_size, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER :: modname = "descr_data"
  INTEGER, PARAMETER :: nproma = 8

  ! public only to ICON
  ! note: some structures defined here and USEd from minimal_example.F90 to fill descriptive data
  !       this circular access pattern is specific to the simplified mockup (combining
  !       the setup of ICON data structures and filling the ComIn structures in one module)
  PUBLIC :: setup_descr_data
  PUBLIC :: update_exposed_descriptive_data
  PUBLIC :: p_patch
  PUBLIC :: nproma, n_dom, n_dom_dt, max_dom, grf_bdywidth_c
  PUBLIC :: min_rlcell_int, min_rlcell, max_rlcell
  PUBLIC :: min_rlvert_int, min_rlvert, max_rlvert
  PUBLIC :: min_rledge_int, min_rledge, max_rledge
  PUBLIC :: grf_bdywidth_e, dt, exp_start, exp_stop, sim_current
  PUBLIC :: run_start, run_stop
  PUBLIC :: expose_descriptive_data
  PUBLIC :: clear_descr_data

  TYPE :: t_grid_cells
    REAL(wp), ALLOCATABLE         :: area(:,:)
    INTEGER                       :: max_connectivity
    INTEGER, ALLOCATABLE          :: num_edges(:,:)
    INTEGER, ALLOCATABLE          :: start_index(:)
    INTEGER, ALLOCATABLE          :: end_index(:)
    INTEGER, ALLOCATABLE          :: start_block(:)
    INTEGER, ALLOCATABLE          :: end_block(:)
    INTEGER, ALLOCATABLE          :: vertex_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: vertex_blk(:,:,:)
    INTEGER, ALLOCATABLE          :: neighbor_blk(:,:,:)
    INTEGER, ALLOCATABLE          :: neighbor_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: edge_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: edge_blk(:,:,:)
    REAL(wp), ALLOCATABLE         :: hhl(:,:,:)
    TYPE(t_geographical_coordinates), ALLOCATABLE :: center(:,:)
    INTEGER, ALLOCATABLE          :: decomp_domain(:,:)
    INTEGER, ALLOCATABLE          :: child_id(:,:)
    INTEGER, ALLOCATABLE          :: child_idx(:,:,:), child_blk(:,:,:)
    INTEGER, ALLOCATABLE          :: parent_glb_idx(:,:), parent_glb_blk(:,:)
  END TYPE t_grid_cells

  TYPE :: t_grid_vertices
    INTEGER, ALLOCATABLE          :: start_index(:)
    INTEGER, ALLOCATABLE          :: end_index(:)
    INTEGER, ALLOCATABLE          :: num_edges(:,:)
    INTEGER, ALLOCATABLE          :: start_block(:)
    INTEGER, ALLOCATABLE          :: end_block(:)
    INTEGER, ALLOCATABLE          :: neighbor_blk(:,:,:)
    INTEGER, ALLOCATABLE          :: neighbor_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: cell_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: cell_blk(:,:,:)
    INTEGER, ALLOCATABLE          :: edge_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: edge_blk(:,:,:)
    TYPE(t_geographical_coordinates), ALLOCATABLE :: vertex(:,:)
  END TYPE t_grid_vertices

  TYPE :: t_grid_edges
    INTEGER, ALLOCATABLE          :: start_index(:)
    INTEGER, ALLOCATABLE          :: end_index(:)
    INTEGER, ALLOCATABLE          :: start_block(:)
    INTEGER, ALLOCATABLE          :: end_block(:)
    INTEGER, ALLOCATABLE          :: cell_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: cell_blk(:,:,:)
    INTEGER, ALLOCATABLE          :: vertex_idx(:,:,:)
    INTEGER, ALLOCATABLE          :: vertex_blk(:,:,:)
    TYPE(t_geographical_coordinates), ALLOCATABLE :: center(:,:)
    INTEGER, ALLOCATABLE          :: child_id(:,:)
    INTEGER, ALLOCATABLE          :: child_idx(:,:,:), child_blk(:,:,:)
    INTEGER, ALLOCATABLE          :: parent_glb_idx(:,:), parent_glb_blk(:,:)
  END TYPE t_grid_edges

  ! from mo_decomposition_tools.f90
  TYPE t_glb2loc_index_lookup
    ! sorted list of global indices.
    ! supposed to be used in a binary_search
    INTEGER, ALLOCATABLE :: glb_index(:)
    ! list of local indices. Supposed to lookup after the global index has beed
    ! found in glb_index
    INTEGER, ALLOCATABLE :: glb_index_to_loc(:)
  END TYPE  t_glb2loc_index_lookup

  TYPE :: t_patch
    CHARACTER(LEN=filename_max)   :: grid_filename
    INTEGER(c_signed_char), ALLOCATABLE        :: grid_uuid(:)
    INTEGER                       :: id
    INTEGER                       :: n_childdom
    INTEGER, ALLOCATABLE          :: child_id(:)
    REAL(wp)                      :: dom_start
    REAL(wp)                      :: dom_end
    INTEGER                       :: nlev
    INTEGER                       :: nshift
    INTEGER                       :: nshift_total
    INTEGER                       :: n_patch_cells
    INTEGER                       :: n_patch_edges
    INTEGER                       :: n_patch_verts
    INTEGER                       :: n_patch_cells_g
    INTEGER                       :: n_patch_edges_g
    INTEGER                       :: n_patch_verts_g
    INTEGER                       :: nblks_c
    INTEGER                       :: nblks_e
    INTEGER                       :: nblks_v
    INTEGER, ALLOCATABLE          :: refin_c_ctrl(:,:)
    INTEGER, ALLOCATABLE          :: refin_e_ctrl(:,:)
    INTEGER, ALLOCATABLE          :: refin_v_ctrl(:,:)
    TYPE(t_grid_cells)            :: cells
    TYPE(t_grid_vertices)         :: verts
    TYPE(t_grid_edges)            :: edges
    INTEGER                       :: global_size
    INTEGER, ALLOCATABLE          :: glb_index(:)
    TYPE (t_glb2loc_index_lookup) :: glb2loc_index
  END TYPE t_patch

  ! not part of p_patch structure
  INTEGER  :: n_dom, n_dom_dt, max_dom
  ! from ./src/shared/mo_impl_constants.f90
  ! INTEGER, PARAMETER :: min_rlcell_int = -4 => see below
  ! from ./src/shared/mo_impl_constants_grf.f90
  INTEGER, PARAMETER :: grf_bdywidth_c = 4
  INTEGER, PARAMETER :: grf_bdywidth_e = 9
  REAL(wp)           :: dt

  CHARACTER(LEN=max_datetime_str_len) :: exp_start
  CHARACTER(LEN=max_datetime_str_len) :: exp_stop
  CHARACTER(LEN=max_datetime_str_len) :: sim_current
  CHARACTER(LEN=max_datetime_str_len) :: run_start
  CHARACTER(LEN=max_datetime_str_len) :: run_stop

  TYPE(t_patch), TARGET, ALLOCATABLE  :: p_patch(:)

  ! from mo_impl_constants
  INTEGER, PARAMETER :: max_hw         = 2                         ! maximum halo width (n_ghost_rows)
  INTEGER, PARAMETER :: min_rlcell_int = -4                        ! previously -6
  INTEGER, PARAMETER :: min_rlcell     = min_rlcell_int - 2*max_hw ! = -8
  INTEGER, PARAMETER :: max_rlcell     = 5                         ! previously 8
  INTEGER, PARAMETER :: min_rlvert_int = min_rlcell_int
  INTEGER, PARAMETER :: min_rlvert     = -4 - (max_hw+1)
  INTEGER, PARAMETER :: max_rlvert     = max_rlcell
  INTEGER, PARAMETER :: min_rledge_int = 2*min_rlcell_int          ! -8
  INTEGER, PARAMETER :: min_rledge     = -8 - (2*max_hw+1)
  INTEGER, PARAMETER :: max_rledge     = 2*max_rlcell

  ! from mo_run_config
  INTEGER, ALLOCATABLE :: number_of_grid_used(:)

  ! module-wide variables for ComIn
  TYPE(t_comin_descrdata_global)                  :: comin_global
  TYPE(t_comin_descrdata_domain), ALLOCATABLE     :: comin_grid(:)
  TYPE(t_comin_descrdata_simulation_interval)     :: comin_simulation_interval

CONTAINS

  !> very rudimentary and incomplete first setup of ICON's descriptive data structures
  SUBROUTINE setup_descr_data(comm)
    INTEGER, INTENT(IN) :: comm
    !
    INTEGER :: status, length, jg, i, mpi_size, mpi_rank, ierr, local_size, jk
    CHARACTER(LEN=128) :: cprefix
    INTEGER, ALLOCATABLE :: blkd_glb_index(:,:)

    CALL MPI_Comm_size(comm, mpi_size, ierr)
    CALL MPI_Comm_rank(comm, mpi_rank, ierr)
    ! set up non p_patch parts
    ! n_dom is the number of domains fully set up in the minimal example
    ! Note: n_dom_dt is the size of a small subpart (n_childdom and child_id)
    ! which has more domains to illustrate the functionality of expose_timesteplength_domain
    n_dom_dt = 5
    n_dom = 1
    max_dom = 10
    dt = 600
    exp_start = "1950-01-01T00:00:00"
    exp_stop = "1950-01-01T00:30:00"
    sim_current = "1950-01-01T00:00:00"
    run_start = "1950-01-01T00:00:00"
    run_stop = "1950-01-01T00:30:00"

    ! set up p_patch
    ALLOCATE(p_patch(n_dom_dt))
    ALLOCATE(number_of_grid_used(n_dom_dt))

    ! read-in of demo grid file
    !
    ! we use the 1st command line argument as a path prefix:
    cprefix = ""
    CALL GET_COMMAND_ARGUMENT (1, cprefix, length, status)

    IF ((status .NE. 0) .OR. (LEN_TRIM(cprefix) == 0)) THEN
      CALL finish("descr_data::setup_descr_data", "Getting command line argument failed!")
    END IF

    ! for now: example setup allocating only cell components
    jg=1
    SELECT CASE(mpi_size)
    CASE(1)
      p_patch(jg)%grid_filename = TRIM(cprefix)//"/grid.1_0.nc"
    CASE(2)
      IF (mpi_rank == 0) THEN
        p_patch(jg)%grid_filename = TRIM(cprefix)//"/grid.2_0.nc"
      ELSE
        p_patch(jg)%grid_filename = TRIM(cprefix)//"/grid.2_1.nc"
      END IF
    CASE DEFAULT
      CALL finish("descr_data::setup_descr_data", "Not implemented!")
    END SELECT

    IF (mpi_rank == 0) THEN
      WRITE (0,*) "     read NetCDF file '", TRIM(p_patch(jg)%grid_filename), "'."
    END IF

    CALL read_netcdf_scalar("clon", nproma, p_patch(jg)%cells%center, "lon",  p_patch(jg)%n_patch_cells, &
         &                             opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("clat", nproma, p_patch(jg)%cells%center, "lat",  &
         &                             opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("vlon", nproma, p_patch(jg)%verts%vertex, "lon", p_patch(jg)%n_patch_verts, &
         &                             opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("vlat", nproma, p_patch(jg)%verts%vertex, "lat", &
         &                             opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("elon", nproma, p_patch(jg)%edges%center, "lon", p_patch(jg)%n_patch_edges, &
         &                             opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("elat", nproma, p_patch(jg)%edges%center, "lat", &
         &                             opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("vertex_of_cell", nproma, p_patch(jg)%cells%vertex_idx, &
         & opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("cells_decomp_domain", nproma, p_patch(jg)%cells%decomp_domain, &
         & opt_filename = p_patch(jg)%grid_filename)
    CALL read_netcdf_scalar("cells_glb_index", 1, blkd_glb_index, &
         & opt_filename = p_patch(jg)%grid_filename)
    p_patch(jg)%glb_index = blkd_glb_index(1,:)

    ! make vertex_idx and vertex_blk like in ICON
    p_patch(jg)%cells%vertex_blk = (p_patch(jg)%cells%vertex_idx -1)/nproma+1
    p_patch(jg)%cells%vertex_idx = MOD(p_patch(jg)%cells%vertex_idx -1, nproma)+1

    local_size = SIZE(p_patch(jg)%glb_index)

    ! compute the global size by summing up owned cells
    p_patch(jg)%global_size = COUNT(p_patch(jg)%cells%decomp_domain .EQ. 0)
    CALL MPI_Allreduce(MPI_IN_PLACE, p_patch(jg)%global_size, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

    ! compute glb2loc lookup by sorting the global indices and storing the perm.
    ALLOCATE(p_patch(jg)%glb2loc_index%glb_index(local_size))
    p_patch(jg)%glb2loc_index%glb_index = p_patch(jg)%glb_index
    ALLOCATE(p_patch(jg)%glb2loc_index%glb_index_to_loc(local_size))
    p_patch(jg)%glb2loc_index%glb_index_to_loc = (/(i, i=1,local_size)/)
    CALL quicksort_permutation_int(p_patch(jg)%glb2loc_index%glb_index, &
                                   p_patch(jg)%glb2loc_index%glb_index_to_loc)

    ! in mock-up only rudimentary setup of p_patch for tests
    DO jg=1,n_dom_dt
      p_patch(jg)%id = jg
      ALLOCATE(p_patch(jg)%child_id(n_dom_dt))
      IF (jg == 1) THEN
        p_patch(jg)%n_childdom = 2
        p_patch(jg)%child_id(1:2) = [2,3]
      ELSE IF (jg==2) THEN
        p_patch(jg)%n_childdom = 1
        p_patch(jg)%child_id(1) = 4
      ELSE IF (jg==3) THEN
        p_patch(jg)%n_childdom = 0
        p_patch(jg)%child_id(1) = 0
      ELSE IF (jg==4) THEN
        p_patch(jg)%n_childdom = 1
        p_patch(jg)%child_id(1) = 5
      ELSE IF (jg==5) THEN
        p_patch(jg)%n_childdom = 0
        p_patch(jg)%child_id(1) = 0
      END IF
    END DO
    DO jg=1, n_dom
      p_patch(jg)%nblks_c = blk_no(p_patch(jg)%n_patch_cells, nproma)
      p_patch(jg)%nblks_e = blk_no(p_patch(jg)%n_patch_edges, nproma)
      p_patch(jg)%nblks_v = blk_no(p_patch(jg)%n_patch_verts, nproma)
      p_patch(jg)%n_patch_cells_g = p_patch(jg)%global_size
      p_patch(jg)%n_patch_edges_g = -666 ! okay. that not true but now its deterministic
      p_patch(jg)%n_patch_verts_g = -666 ! okay. that not true but now its deterministic
      p_patch(jg)%nlev = 10
      p_patch(jg)%dom_start = 0.0
      p_patch(jg)%dom_end = 1.0
      p_patch(jg)%nshift = 42
      p_patch(jg)%nshift_total = 43
      p_patch(jg)%grid_uuid = [42_1]
      ALLOCATE(p_patch(jg)%cells%start_index(min_rlcell:max_rlcell))
      p_patch(jg)%cells%start_index(1) = 4
      ALLOCATE(p_patch(jg)%cells%end_index(min_rlcell:max_rlcell))
      p_patch(jg)%cells%end_index(1) = 3
      ALLOCATE(p_patch(jg)%cells%start_block(min_rlcell:max_rlcell))
      ALLOCATE(p_patch(jg)%cells%end_block(min_rlcell:max_rlcell))
      ALLOCATE(p_patch(jg)%cells%neighbor_blk(nproma,p_patch(jg)%nblks_c,3))
      ALLOCATE(p_patch(jg)%cells%neighbor_idx(nproma,p_patch(jg)%nblks_c,3))
      ALLOCATE(p_patch(jg)%cells%edge_idx(nproma,p_patch(jg)%nblks_c,3))
      ALLOCATE(p_patch(jg)%cells%edge_blk(nproma,p_patch(jg)%nblks_c,3))
      ALLOCATE(p_patch(jg)%cells%area(nproma,p_patch(jg)%nblks_c))
      ALLOCATE(p_patch(jg)%cells%hhl(nproma,p_patch(jg)%nlev+1,p_patch(jg)%nblks_c))
      p_patch(jg)%cells%hhl(:,p_patch(jg)%nlev+1,:) = 0._wp
      DO jk = p_patch(jg)%nlev,1,-1
        p_patch(jg)%cells%hhl(:,jk,:) = &
          p_patch(jg)%cells%hhl(:,jk+1,:) + 100._wp
      END DO

      ALLOCATE(p_patch(jg)%cells%num_edges(nproma,p_patch(jg)%nblks_c))
      ALLOCATE(p_patch(jg)%cells%child_id(nproma,p_patch(jg)%nblks_c))
      ALLOCATE(p_patch(jg)%cells%child_idx(nproma,p_patch(jg)%nblks_c,4))
      ALLOCATE(p_patch(jg)%cells%child_blk(nproma,p_patch(jg)%nblks_c,4))
      ALLOCATE(p_patch(jg)%cells%parent_glb_idx(nproma,p_patch(jg)%nblks_c))
      ALLOCATE(p_patch(jg)%cells%parent_glb_blk(nproma,p_patch(jg)%nblks_c))
      number_of_grid_used(jg) = 2
      ALLOCATE(p_patch(jg)%refin_c_ctrl(nproma, p_patch(jg)%nblks_c))
      ALLOCATE(p_patch(jg)%refin_e_ctrl(nproma, p_patch(jg)%nblks_c))
      ALLOCATE(p_patch(jg)%refin_v_ctrl(nproma, p_patch(jg)%nblks_c))

      ALLOCATE(p_patch(jg)%edges%start_index(min_rledge:max_rledge))
      p_patch(jg)%edges%start_index(1) = 4
      ALLOCATE(p_patch(jg)%edges%end_index(min_rledge:max_rledge))
      p_patch(jg)%edges%end_index(1) = 3
      ALLOCATE(p_patch(jg)%edges%start_block(min_rledge:max_rledge))
      ALLOCATE(p_patch(jg)%edges%end_block(min_rledge:max_rledge))
      ALLOCATE(p_patch(jg)%edges%child_id(nproma,p_patch(jg)%nblks_e))
      ALLOCATE(p_patch(jg)%edges%child_idx(nproma,p_patch(jg)%nblks_e,4))
      ALLOCATE(p_patch(jg)%edges%child_blk(nproma,p_patch(jg)%nblks_e,4))
      ALLOCATE(p_patch(jg)%edges%parent_glb_idx(nproma,p_patch(jg)%nblks_e))
      ALLOCATE(p_patch(jg)%edges%parent_glb_blk(nproma,p_patch(jg)%nblks_e))
      ALLOCATE(p_patch(jg)%edges%cell_idx(nproma,p_patch(jg)%nblks_e,2))
      ALLOCATE(p_patch(jg)%edges%cell_blk(nproma,p_patch(jg)%nblks_e,2))
      ALLOCATE(p_patch(jg)%edges%vertex_idx(nproma,p_patch(jg)%nblks_e,2))
      ALLOCATE(p_patch(jg)%edges%vertex_blk(nproma,p_patch(jg)%nblks_e,2))
      ALLOCATE(p_patch(jg)%verts%num_edges(nproma,p_patch(jg)%nblks_v))
      ALLOCATE(p_patch(jg)%verts%start_index(min_rlvert:max_rlvert))
      p_patch(jg)%verts%start_index(1) = 4
      ALLOCATE(p_patch(jg)%verts%end_index(min_rlvert:max_rlvert))
      p_patch(jg)%verts%end_index(1) = 3
      ALLOCATE(p_patch(jg)%verts%start_block(min_rlvert:max_rlvert))
      ALLOCATE(p_patch(jg)%verts%end_block(min_rlvert:max_rlvert))
      ALLOCATE(p_patch(jg)%verts%cell_idx(nproma,p_patch(jg)%nblks_v,6))
      ALLOCATE(p_patch(jg)%verts%cell_blk(nproma,p_patch(jg)%nblks_v,6))
      ALLOCATE(p_patch(jg)%verts%edge_idx(nproma,p_patch(jg)%nblks_c,6))
      ALLOCATE(p_patch(jg)%verts%edge_blk(nproma,p_patch(jg)%nblks_c,6))
      ALLOCATE(p_patch(jg)%verts%neighbor_blk(nproma,p_patch(jg)%nblks_c,6))
      ALLOCATE(p_patch(jg)%verts%neighbor_idx(nproma,p_patch(jg)%nblks_c,6))
    END DO
  END SUBROUTINE setup_descr_data

  !> General wrapper subroutine for exposing descriptive data
  !> structures to ComIn.
  SUBROUTINE expose_descriptive_data(n_dom, max_dom, nproma, min_rlcell_int,   &
       &                  min_rlcell, max_rlcell, min_rlvert_int, min_rlvert,  &
       &                  max_rlvert, min_rledge_int, min_rledge, max_rledge,  &
       &                  grf_bdywidth_c, grf_bdywidth_e,          &
       &                  dt, lrestart, sim_time_start, sim_time_end,          &
       &                  sim_time_current, run_time_start,                    &
       &                  run_time_stop, yac_instance_id)
    INTEGER,                   INTENT(IN) :: n_dom
    INTEGER,                   INTENT(IN) :: max_dom
    INTEGER,                   INTENT(IN) :: nproma
    INTEGER,                   INTENT(IN) :: min_rlcell_int
    INTEGER,                   INTENT(IN) :: min_rlcell
    INTEGER,                   INTENT(IN) :: max_rlcell
    INTEGER,                   INTENT(IN) :: min_rlvert_int
    INTEGER,                   INTENT(IN) :: min_rlvert
    INTEGER,                   INTENT(IN) :: max_rlvert
    INTEGER,                   INTENT(IN) :: min_rledge_int
    INTEGER,                   INTENT(IN) :: min_rledge
    INTEGER,                   INTENT(IN) :: max_rledge
    INTEGER,                   INTENT(IN) :: grf_bdywidth_c
    INTEGER,                   INTENT(IN) :: grf_bdywidth_e
    REAL(wp),                  INTENT(IN) :: dt
    LOGICAL,                   INTENT(IN) :: lrestart
    CHARACTER(LEN=*),          INTENT(IN) :: sim_time_start
    CHARACTER(LEN=*),          INTENT(IN) :: sim_time_end
    CHARACTER(LEN=*),          INTENT(IN) :: sim_time_current
    CHARACTER(LEN=*),          INTENT(IN) :: run_time_start
    CHARACTER(LEN=*),          INTENT(IN) :: run_time_stop
    INTEGER,                   INTENT(IN) :: yac_instance_id
    !

    ! global data need to be transferred first - required for grid
    ! information
    CALL expose_global(n_dom, max_dom, nproma, min_rlcell_int, min_rlcell, &
         &             max_rlcell, min_rlvert_int, min_rlvert, max_rlvert, &
         &             min_rledge_int, min_rledge, max_rledge,             &
         &             grf_bdywidth_c, grf_bdywidth_e, lrestart, yac_instance_id)
    CALL expose_grid_info(p_patch)
    CALL expose_simulation_interval(sim_time_start, sim_time_end, sim_time_current, &
             &             run_time_start, run_time_stop)
    CALL comin_descrdata_set_fct_glb2loc_cell(glb2loc_index)
    CALL expose_timesteplength_domain(1, dt)
  END SUBROUTINE expose_descriptive_data

  !> Update of exposed descriptive data (ComIn)
  SUBROUTINE update_exposed_descriptive_data(sim_time_current)
    CHARACTER(LEN=*),          INTENT(IN) :: sim_time_current

    CALL comin_current_set_datetime(sim_time_current)
  END SUBROUTINE update_exposed_descriptive_data

  !> fill ComIn descriptive data structures from ICON
  SUBROUTINE expose_global(n_dom, max_dom, nproma, min_rlcell_int, min_rlcell, &
       &                   max_rlcell, min_rlvert_int, min_rlvert, max_rlvert, &
       &                   min_rledge_int, min_rledge, max_rledge,             &
       &                   grf_bdywidth_c, grf_bdywidth_e, lrestart, yac_instance_id)
    INTEGER, INTENT(IN) :: n_dom
    INTEGER, INTENT(IN) :: max_dom
    INTEGER, INTENT(IN) :: nproma
    INTEGER, INTENT(IN) :: min_rlcell_int
    INTEGER, INTENT(IN) :: min_rlcell
    INTEGER, INTENT(IN) :: max_rlcell
    INTEGER, INTENT(IN) :: min_rlvert_int
    INTEGER, INTENT(IN) :: min_rlvert
    INTEGER, INTENT(IN) :: max_rlvert
    INTEGER, INTENT(IN) :: min_rledge_int
    INTEGER, INTENT(IN) :: min_rledge
    INTEGER, INTENT(IN) :: max_rledge
    INTEGER, INTENT(IN) :: grf_bdywidth_c
    INTEGER, INTENT(IN) :: grf_bdywidth_e
    LOGICAL, INTENT(IN) :: lrestart
    INTEGER, INTENT(IN) :: yac_instance_id

    comin_global%n_dom            = n_dom
    comin_global%max_dom          = max_dom
    comin_global%nproma           = nproma
    comin_global%min_rlcell_int   = min_rlcell_int
    comin_global%min_rlcell       = min_rlcell
    comin_global%max_rlcell       = max_rlcell
    comin_global%min_rlvert_int   = min_rlvert_int
    comin_global%min_rlvert       = min_rlvert
    comin_global%max_rlvert       = max_rlvert
    comin_global%min_rledge_int   = min_rledge_int
    comin_global%min_rledge       = min_rledge
    comin_global%max_rledge       = max_rledge
    comin_global%grf_bdywidth_c   = grf_bdywidth_c
    comin_global%grf_bdywidth_e   = grf_bdywidth_e
    comin_global%lrestartrun      = lrestart
    comin_global%vct_a            = [42.] ! what else?!
    comin_global%yac_instance_id  = yac_instance_id
    comin_global%host_git_remote_url = "dummy_url"
    comin_global%host_git_branch  = "dummy_branch"
    comin_global%host_git_tag     = "dummy_tag"
    comin_global%host_revision    = "dummy_revision"
    comin_global%device_name      = "Test device"
    comin_global%device_vendor    = "Dummy vendor"
    comin_global%device_driver    = "Dummy driver"

    ! register global info in ComIn
    CALL comin_descrdata_set_global(comin_global)
  END SUBROUTINE expose_global

  !> fill ComIn descriptive data structures from ICON
  SUBROUTINE expose_grid_info(patch)
    TYPE(t_patch), TARGET,     INTENT(IN)      :: patch(:)
    !local var
    INTEGER :: jg

    ! note that this prototype implementation is missing various members of the ComIn
    ! descriptive data structures
    ALLOCATE(comin_grid(SIZE(patch)))
    DO jg=1,n_dom
      ALLOCATE(comin_grid(jg)%number_of_grid_used(SIZE(patch)))
      comin_grid(jg)%dom_start = patch(jg)%dom_start
      comin_grid(jg)%dom_end = patch(jg)%dom_end
      ! cell, vertex and edge independent variables
      comin_grid(jg)%grid_filename => patch(jg)%grid_filename
      comin_grid(jg)%grid_uuid =>  patch(jg)%grid_uuid
      comin_grid(jg)%number_of_grid_used = number_of_grid_used
      comin_grid(jg)%id => patch(jg)%id
      comin_grid(jg)%n_childdom => patch(jg)%n_childdom
      comin_grid(jg)%nlev => patch(jg)%nlev
      comin_grid(jg)%nshift => patch(jg)%nshift
      comin_grid(jg)%nshift_total => patch(jg)%nshift_total
      comin_grid(jg)%cells%ncells => patch(jg)%n_patch_cells
      comin_grid(jg)%cells%ncells_global => patch(jg)%n_patch_cells_g
      comin_grid(jg)%edges%nedges => patch(jg)%n_patch_edges
      comin_grid(jg)%edges%nedges_global => patch(jg)%n_patch_edges_g
      comin_grid(jg)%verts%nverts => patch(jg)%n_patch_verts
      comin_grid(jg)%verts%nverts_global => patch(jg)%n_patch_verts_g
      comin_grid(jg)%cells%nblks => patch(jg)%nblks_c
      comin_grid(jg)%edges%nblks => patch(jg)%nblks_e
      comin_grid(jg)%verts%nblks => patch(jg)%nblks_v
      comin_grid(jg)%cells%refin_ctrl =>  patch(jg)%refin_c_ctrl
      comin_grid(jg)%edges%refin_ctrl =>  patch(jg)%refin_e_ctrl
      comin_grid(jg)%verts%refin_ctrl =>  patch(jg)%refin_v_ctrl
      ! cell components
      comin_grid(jg)%cells%max_connectivity => patch(jg)%cells%max_connectivity
      comin_grid(jg)%cells%num_edges => patch(jg)%cells%num_edges
      comin_grid(jg)%cells%start_index => patch(jg)%cells%start_index
      comin_grid(jg)%cells%end_index => patch(jg)%cells%end_index
      comin_grid(jg)%cells%start_block => patch(jg)%cells%start_block
      comin_grid(jg)%cells%end_block => patch(jg)%cells%end_block
      comin_grid(jg)%cells%vertex_idx => patch(jg)%cells%vertex_idx
      comin_grid(jg)%cells%vertex_blk => patch(jg)%cells%vertex_blk
      comin_grid(jg)%cells%neighbor_idx => patch(jg)%cells%neighbor_idx
      comin_grid(jg)%cells%neighbor_blk => patch(jg)%cells%neighbor_blk
      comin_grid(jg)%cells%edge_idx => patch(jg)%cells%edge_idx
      comin_grid(jg)%cells%edge_blk => patch(jg)%cells%edge_blk
      comin_grid(jg)%cells%clon = patch(jg)%cells%center%lon
      comin_grid(jg)%cells%clat = patch(jg)%cells%center%lat
      comin_grid(jg)%cells%area => patch(jg)%cells%area
      comin_grid(jg)%cells%hhl = patch(jg)%cells%hhl
      comin_grid(jg)%cells%glb_index => patch(jg)%glb_index
      comin_grid(jg)%cells%decomp_domain => patch(jg)%cells%decomp_domain
      comin_grid(jg)%cells%child_id  => patch(jg)%cells%child_id
      comin_grid(jg)%cells%child_idx => patch(jg)%cells%child_idx
      comin_grid(jg)%cells%child_blk => patch(jg)%cells%child_blk
      comin_grid(jg)%cells%parent_glb_idx => patch(jg)%cells%parent_glb_idx
      comin_grid(jg)%cells%parent_glb_blk => patch(jg)%cells%parent_glb_blk
      ! edge components
      comin_grid(jg)%edges%child_id => patch(jg)%edges%child_id
      comin_grid(jg)%edges%child_idx => patch(jg)%edges%child_idx
      comin_grid(jg)%edges%child_blk => patch(jg)%edges%child_blk
      comin_grid(jg)%edges%cell_idx => patch(jg)%edges%cell_idx
      comin_grid(jg)%edges%cell_blk => patch(jg)%edges%cell_blk
      comin_grid(jg)%edges%elon = patch(jg)%edges%center%lon
      comin_grid(jg)%edges%elat = patch(jg)%edges%center%lat
      comin_grid(jg)%edges%start_block => patch(jg)%edges%start_block
      comin_grid(jg)%edges%end_block => patch(jg)%edges%end_block
      comin_grid(jg)%edges%start_index => patch(jg)%edges%start_index
      comin_grid(jg)%edges%end_index => patch(jg)%edges%end_index
      comin_grid(jg)%edges%parent_glb_idx => patch(jg)%edges%parent_glb_idx
      comin_grid(jg)%edges%parent_glb_blk => patch(jg)%edges%parent_glb_blk
      comin_grid(jg)%edges%vertex_idx => patch(jg)%edges%vertex_idx
      comin_grid(jg)%edges%vertex_blk => patch(jg)%edges%vertex_blk
      ! vertex components
      comin_grid(jg)%verts%vlon = patch(jg)%verts%vertex%lon
      comin_grid(jg)%verts%vlat = patch(jg)%verts%vertex%lat
      comin_grid(jg)%verts%num_edges => patch(jg)%verts%num_edges
      comin_grid(jg)%verts%cell_idx => patch(jg)%verts%cell_idx
      comin_grid(jg)%verts%cell_blk => patch(jg)%verts%cell_blk
      comin_grid(jg)%verts%start_block => patch(jg)%verts%start_block
      comin_grid(jg)%verts%end_block => patch(jg)%verts%end_block
      comin_grid(jg)%verts%start_index => patch(jg)%verts%start_index
      comin_grid(jg)%verts%end_index => patch(jg)%verts%end_index
      comin_grid(jg)%verts%edge_idx => patch(jg)%verts%edge_idx
      comin_grid(jg)%verts%edge_blk => patch(jg)%verts%edge_blk
      comin_grid(jg)%verts%neighbor_idx => patch(jg)%verts%neighbor_idx
      comin_grid(jg)%verts%neighbor_blk => patch(jg)%verts%neighbor_blk
    END DO

    ! register grid info in ComIn
    CALL comin_descrdata_set_domain(comin_grid)
  END SUBROUTINE expose_grid_info

  SUBROUTINE expose_simulation_interval(exp_start, exp_stop, sim_current, run_start, run_stop)
    CHARACTER(LEN=*), INTENT(IN) :: exp_start
    CHARACTER(LEN=*), INTENT(IN) :: exp_stop
    CHARACTER(LEN=*), INTENT(IN) :: sim_current
    CHARACTER(LEN=*), INTENT(IN) :: run_start
    CHARACTER(LEN=*), INTENT(IN) :: run_stop

    comin_simulation_interval%exp_start   = exp_start
    comin_simulation_interval%exp_stop     = exp_stop
    comin_simulation_interval%run_start   = run_start
    comin_simulation_interval%run_stop    = run_stop

    ! register time info in ComIn
    CALL comin_descrdata_set_simulation_interval(comin_simulation_interval)
    CALL comin_current_set_datetime(sim_current)

  END SUBROUTINE expose_simulation_interval

  RECURSIVE SUBROUTINE expose_timesteplength_domain(jg, dt_current)
    INTEGER,  INTENT(IN) :: jg
    REAL(wp), INTENT(IN) :: dt_current
    !LOCAL
    INTEGER :: jn

    CALL comin_descrdata_set_timesteplength(jg, dt_current)

    DO jn=1,p_patch(jg)%n_childdom
      CALL expose_timesteplength_domain(p_patch(jg)%child_id(jn), dt_current/2.0_wp)
    END DO
  END SUBROUTINE expose_timesteplength_domain

  SUBROUTINE clear_descr_data()
    ! do nothing
  END SUBROUTINE clear_descr_data

  ! taken from ICONs `mo_util_sort.f90`:
  SUBROUTINE swap_permutation_int(a, i,j, permutation)
    !> array for in-situ sorting
    INTEGER,  INTENT(INOUT)           :: a(:)
    !> indices to be exchanged
    INTEGER,  INTENT(IN)              :: i,j
    !> (optional) permutation of indices
    INTEGER,  INTENT(INOUT)           :: permutation(:)
    ! local variables
    INTEGER :: t, t_p

    t    = a(i)
    a(i) = a(j)
    a(j) = t
    t_p            = permutation(i)
    permutation(i) = permutation(j)
    permutation(j) = t_p
  END SUBROUTINE swap_permutation_int

  ! --------------------------------------------------------------------
  !> Simple recursive implementation of Hoare's QuickSort algorithm
  !  for a 1D array of INTEGER values.
  !
  !  Ordering after the sorting process: smallest...largest.
  !
  RECURSIVE SUBROUTINE quicksort_permutation_int(a, permutation, l_in, r_in)
    INTEGER,  INTENT(INOUT)           :: a(:)           !< array for in-situ sorting
    INTEGER,  INTENT(INOUT)           :: permutation(:) !< (optional) permutation of indices
    INTEGER,  INTENT(IN),    OPTIONAL :: l_in,r_in      !< left, right partition indices
    ! local variables
    INTEGER :: i,j,l,r,t_p,t,v,m

    IF (PRESENT(l_in)) THEN
      l = l_in
    ELSE
      l = 1
    END IF
    IF (PRESENT(r_in)) THEN
      r = r_in
    ELSE
      r = SIZE(a,1)
    END IF
    IF (r>l) THEN
      i = l-1
      j = r

      ! median-of-three selection of partitioning element
      IF ((r-l) > 3) THEN
        m = (l+r)/2
        IF (a(l)>a(m))  CALL swap_permutation_int(a, l,m, permutation)
        IF (a(l)>a(r)) THEN
          CALL swap_permutation_int(a, l,r, permutation)
        ELSE IF (a(r)>a(m)) THEN
          CALL swap_permutation_int(a, r,m, permutation)
        END IF
      END IF

      v = a(r)
      LOOP : DO
        CNTLOOP1 : DO
          i = i+1
          IF (a(i) >= v) EXIT CNTLOOP1
        END DO CNTLOOP1
        CNTLOOP2 : DO
          j = j-1
          IF ((a(j) <= v) .OR. (j==1)) EXIT CNTLOOP2
        END DO CNTLOOP2
        t    = a(i)
        a(i) = a(j)
        a(j) = t

        t_p            = permutation(i)
        permutation(i) = permutation(j)
        permutation(j) = t_p

        IF (j <= i) EXIT LOOP
      END DO LOOP
      a(j) = a(i)
      a(i) = a(r)
      a(r) = t

      permutation(j) = permutation(i)
      permutation(i) = permutation(r)
      permutation(r) = t_p
      CALL quicksort_permutation_int(a,permutation,l,i-1)
      CALL quicksort_permutation_int(a,permutation,i+1,r)
    END IF
  END SUBROUTINE quicksort_permutation_int

  ! taken from ICONs mo_decomposition_tools.f90
  PURE FUNCTION binary_search(array, key)

    INTEGER, INTENT(IN) :: array(:), key
    INTEGER :: binary_search

    INTEGER :: lb, ub, middle

    !$ACC ROUTINE SEQ

    lb = 1
    ub = SIZE(array)
    middle = ub / 2

    IF (ub == 0) THEN
      binary_search = 0
      RETURN
    END IF

    DO WHILE (ub >= lb)

      middle = (ub + lb) / 2
      IF (array(middle) < key) THEN
        lb = middle + 1
      ELSE IF (array(middle) > key) THEN
        ub = middle - 1
      ELSE
        EXIT
      END IF
    END DO

    IF (array(middle) == key) THEN
      binary_search = middle
    ELSE IF (array(middle) > key) THEN
      binary_search = -middle + 1
    ELSE
      binary_search = -middle
    END IF
  END FUNCTION binary_search

  INTEGER FUNCTION glb2loc_index(jg, glb) RESULT(loc)
    INTEGER, INTENT(IN) :: jg
    INTEGER, INTENT(IN) :: glb
    INTEGER :: pos
    pos = binary_search(p_patch(jg)%glb2loc_index%glb_index, glb)
    loc = p_patch(jg)%glb2loc_index%glb_index_to_loc(pos)
  END FUNCTION glb2loc_index

END MODULE descr_data
