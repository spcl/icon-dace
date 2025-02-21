! Contains the definition of coefficients used for div-grad-curl and reconstruction/scalar product.
! All coefficients are three-dimensional arrays including the number of vertical levels. This is necessary
! if one has coefficients that vary within the vertical level but not in time such that one can precompute the
! coefficients. This is in the ocean model where the land-sea mask is different at each level and therefore
! the expansion coefficients associated with land and boundary vary with the vertical level but are constant in time.
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

!----------------------------
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

!----------------------------
MODULE mo_operator_ocean_coeff_3d
  !-------------------------------------------------------------------------

  USE mo_kind,                ONLY: wp, sp
  USE mo_impl_constants,      ONLY: success,&
    &                               max_char_length, beta_plane_coriolis,full_coriolis, &
    &                               SEA_BOUNDARY, BOUNDARY, SEA, min_dolic
  USE mo_math_constants,      ONLY: deg2rad, pi!, rad2deg
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_math_utilities,      ONLY: gc2cc, cc2gc, vector_product, &
    &  arc_length, cvec2gvec
  USE mo_ocean_nml,           ONLY: n_zlev, no_tracer, &
    & coriolis_type, basin_center_lat, basin_height_deg, &
    & select_solver, select_gmres_mp_r, &
    & select_lhs, select_lhs_matrix
  USE mo_exception,           ONLY: message, finish
  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_parallel_config,     ONLY: nproma
  USE mo_sync,                ONLY: sync_c, sync_e, sync_v, sync_patch_array!, sync_idx, global_max
  USE mo_ocean_types,         ONLY: t_hydro_ocean_state, t_operator_coeff, &
    & t_solverCoeff_singlePrecision, t_verticalAdvection_ppm_coefficients
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range
  USE mo_grid_config,         ONLY: grid_sphere_radius, grid_angular_velocity
  USE mo_run_config,          ONLY: dtime
  USE mo_cf_convention,       ONLY: t_cf_var
  USE mo_grib2,               ONLY: t_grib2_var
  USE mo_util_dbg_prnt,       ONLY: dbg_print
  USE mo_grid_geometry_info
  USE mo_fortran_tools,       ONLY: set_acc_host_or_device

  IMPLICIT NONE


  PRIVATE

  PUBLIC  :: t_operator_coeff
  PUBLIC  :: construct_operators_coefficients
  PUBLIC  :: destruct_operators_coefficients
  PUBLIC  :: Get3DVectorToPlanarLocal, Get3DVectorTo2DLocal_array3D
  ! PUBLIC  :: update_diffusion_matrices


  !these two parameters are set below in sbr "allocate_operators_coefficients"
  !according to MAXVAL(patch_2D%cells%num_edges) and MAXVAL(patch_2D%verts%num_edges)
  INTEGER,PUBLIC :: no_dual_edges
  INTEGER,PUBLIC :: no_primal_edges 

  ! flags for computing ocean coefficients
  LOGICAL, PARAMETER :: MID_POINT_DUAL_EDGE = .TRUE. !Please do not change this unless you are sure, you know what you do.
  LOGICAL, PARAMETER :: LARC_LENGTH = .FALSE.
  CHARACTER(LEN=*), PARAMETER :: this_mod_name = 'opcoeff'
  INTEGER :: idt_src    = 1               ! Level of detail for 1 line debug


    REAL(wp), POINTER   :: prime_edge_length      (:,:)
    REAL(wp), POINTER   :: dual_edge_length       (:,:)
    REAL(wp), POINTER   :: dist_cell2edge         (:,:,:)
    TYPE(t_cartesian_coordinates), POINTER :: dual_edge_middle(:,:)

CONTAINS
  !-------------------------------------------------------------------------
  SUBROUTINE init_geometry_forCoefficients(patch_2D)
    TYPE(t_patch), TARGET, INTENT(INOUT)     :: patch_2D

    TYPE(t_subset_range), POINTER :: owned_edges         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_cells         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_verts         ! these are the owned entities
    REAL(wp) :: inverse_sphere_radius
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
    INTEGER :: vertex_1_index, vertex_1_block, vertex_2_index, vertex_2_block
    INTEGER :: return_status
    !-----------------------------------------------------------------------
    inverse_sphere_radius = 1.0_wp / grid_sphere_radius

    owned_edges => patch_2D%edges%owned
    owned_cells => patch_2D%cells%owned
    owned_verts => patch_2D%verts%owned

    ALLOCATE(prime_edge_length(1:nproma,patch_2D%nblks_e), &
      &      dual_edge_length (1:nproma,patch_2D%nblks_e), &
      &      dist_cell2edge   (1:nproma,1:patch_2D%nblks_e,1:2), &
      & stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('init_geometry_forCoefficients', 'allocation failed')
    ENDIF

    ! compute some basic distances
    ! this is required if the cartesian distance is used
    ! instead of the spherical
    !
    ! computes_dist_cell2edge( patch_2D, intp_2D_coeff)
    !
    IF ( MID_POINT_DUAL_EDGE ) THEN
      dual_edge_middle => patch_2D%edges%cartesian_dual_middle
    ELSE
      dual_edge_middle => patch_2D%edges%cartesian_center
    ENDIF

    ! 1) calcultate prima and dual length as cartesian distance
    IF (LARC_LENGTH) THEN

      ! 1a) we just need to get them from the grid
      ! NOTE:  these are earth's distances, translate on a unit sphere
      dist_cell2edge(:,:,:) = &
        & patch_2D%edges%edge_cell_length(:,:,:) * inverse_sphere_radius
      prime_edge_length(:,:) = &
        & patch_2D%edges%primal_edge_length(:,:) * inverse_sphere_radius
      dual_edge_length(:,:) = &
        & patch_2D%edges%dual_edge_length(:,:) * inverse_sphere_radius

    ELSE

      !1b) calcultate prima and dual length as cartesian distance
      prime_edge_length(:,:) = 0.0_wp
      dual_edge_length (:,:) = 0.0_wp
      dist_cell2edge (:,:,:) = 0.0_wp

      DO edge_block = owned_edges%start_block, owned_edges%end_block
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index = start_index, end_index

          !----------------------------------------
          ! calculate the cartesian edge length
          vertex_1_index = patch_2D%edges%vertex_idx(edge_index, edge_block, 1)
          vertex_1_block = patch_2D%edges%vertex_blk(edge_index, edge_block, 1)
          vertex_2_index = patch_2D%edges%vertex_idx(edge_index, edge_block, 2)
          vertex_2_block = patch_2D%edges%vertex_blk(edge_index, edge_block, 2)

          prime_edge_length(edge_index,edge_block) = planar_distance( &
            & patch_2D%verts%cartesian(vertex_1_index, vertex_1_block), &
            & patch_2D%verts%cartesian(vertex_2_index, vertex_2_block), &
            & patch_2D%geometry_info)
          !----------------------------------------

          !----------------------------------------
          ! calculate the cartesian distance of the edge center to the cell center
          DO neigbor = 1,2

            dist_cell2edge(edge_index,edge_block,neigbor) = 0.0_wp

            cell_index = patch_2D%edges%cell_idx(edge_index,edge_block,neigbor)
            cell_block = patch_2D%edges%cell_blk(edge_index,edge_block,neigbor)

            IF (cell_index > 0) THEN
              dist_cell2edge(edge_index,edge_block,neigbor) = planar_distance( &
                & patch_2D%edges%cartesian_center(edge_index,edge_block),    &
                & patch_2D%cells%cartesian_center(cell_index,cell_block),    &
                & patch_2D%geometry_info)
            ENDIF

          ENDDO ! neigbor = 1,2
          !----------------------------------------

          !----------------------------------------
          ! calculate the cartesian dual edge length
          cell_1_index = patch_2D%edges%cell_idx(edge_index, edge_block, 1)
          cell_1_block = patch_2D%edges%cell_blk(edge_index, edge_block, 1)
          cell_2_index = patch_2D%edges%cell_idx(edge_index, edge_block, 2)
          cell_2_block = patch_2D%edges%cell_blk(edge_index, edge_block, 2)

!           IF (cell_1_index > 0 .AND. cell_2_index > 0) THEN
!
!             dual_edge_length(edge_index,edge_block) = planar_distance(         &
!               & patch_2D%cells%cartesian_center(cell_1_index, cell_1_block), &
!               & patch_2D%cells%cartesian_center(cell_2_index, cell_2_block), &
!               & patch_2D%geometry_info)
!
!           ELSE
              dual_edge_length(edge_index,edge_block) =              &
                & dist_cell2edge(edge_index,edge_block,1) + &
                & dist_cell2edge(edge_index,edge_block,2)
!            ENDIF
          !----------------------------------------

        ENDDO ! edge_index=start_index,end_index
      ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

      ! synchronize the edge distances
      CALL sync_patch_array(SYNC_E, patch_2D, dist_cell2edge(:,:,1))
      CALL sync_patch_array(SYNC_E, patch_2D, dist_cell2edge(:,:,2))
      CALL sync_patch_array(SYNC_E, patch_2D, prime_edge_length(:,:))
      CALL sync_patch_array(SYNC_E, patch_2D, dual_edge_length(:,:))
    ENDIF
    ! primal end dual edge lenght have been computed
    !-------------------------------------------

  END SUBROUTINE init_geometry_forCoefficients
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE destruct_geometry_forCoefficients()

    DEALLOCATE(prime_edge_length, &
      &        dual_edge_length,    &
      &        dist_cell2edge)

  END SUBROUTINE destruct_geometry_forCoefficients
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE Get3DVectorTo2DLocal_array3D(vector, position_local, levels, subset, geometry_info, x, y, lacc)
    TYPE(t_cartesian_coordinates), POINTER :: vector(:,:,:)
    TYPE(t_geographical_coordinates) , TARGET :: position_local(:,:)
    INTEGER, POINTER :: levels(:,:) 
    TYPE(t_subset_range), POINTER :: subset 
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info    
    REAL(wp), POINTER ::  x(:,:,:), y(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: lacc
    
    INTEGER  :: blockNo, start_index, end_index, this_index, level
    REAL(wp) :: sinLon, cosLon, sinLat, cosLat
    REAL(wp) :: cartesian_x, cartesian_y, cartesian_z, y_help
    LOGICAL  :: lzacc
    
    CHARACTER(LEN=*), PARAMETER :: method_name='Get3DVectorTo2DLocal_array3D'

    CALL set_acc_host_or_device(lzacc, lacc)
    
    SELECT CASE(geometry_info%geometry_type)

!     CASE (planar_torus_geometry)
!       CALL finish(method_name, "planar_torus_geometry is not implemented yet")
      
    CASE (sphere_geometry)
    
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, this_index, level, sinLon, cosLon, &
!ICON_OMP sinLat, cosLat, cartesian_x, cartesian_y, cartesian_z, y_help) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = subset%start_block, subset%end_block
        CALL get_index_range(subset, blockNo, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO this_index =  start_index, end_index
        
          ! these should be calclulated once and stored in the coefficients structure
          sinLon = SIN(position_local(this_index,blockNo)%lon)
          cosLon = COS(position_local(this_index,blockNo)%lon)
          sinLat = SIN(position_local(this_index,blockNo)%lat)
          cosLat = COS(position_local(this_index,blockNo)%lat)

          !$ACC LOOP SEQ
          DO level = 1, levels(this_index,blockNo)
            cartesian_x = vector(this_index,level,blockNo)%x(1)
            cartesian_y = vector(this_index,level,blockNo)%x(2)
            cartesian_z = vector(this_index,level,blockNo)%x(3)
            
            x(this_index,level,blockNo) = cosLon * cartesian_y - sinLon * cartesian_x
            y_help = cosLon * cartesian_x + sinLon * cartesian_y
            y_help = sinLat * y_help
            y(this_index,level,blockNo) = cosLat * cartesian_z - y_help
                        
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
            
              
    CASE ( planar_channel_geometry,  planar_geometry, planar_torus_geometry)
    
      ! just a projection
!ICON_OMP_PARALLEL_DO PRIVATE(start_index,end_index, this_index, level) ICON_OMP_DEFAULT_SCHEDULE
      DO blockNo = subset%start_block, subset%end_block
        CALL get_index_range(subset, blockNo, start_index, end_index)
        !$ACC PARALLEL LOOP GANG VECTOR DEFAULT(PRESENT) ASYNC(1) IF(lzacc)
        DO this_index =  start_index, end_index
          DO level = 1, levels(this_index,blockNo)  
            x(this_index,level,blockNo) = vector(this_index,level,blockNo)%x(1)
            y(this_index,level,blockNo) = vector(this_index,level,blockNo)%x(2)
          ENDDO
        ENDDO
        !$ACC END PARALLEL LOOP
      ENDDO
      !$ACC WAIT(1)
!ICON_OMP_END_PARALLEL_DO
      
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT
  END SUBROUTINE Get3DVectorTo2DLocal_array3D
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE Get3DVectorToPlanarLocal(vector, position_local, geometry_info, x, y)
    TYPE(t_cartesian_coordinates), INTENT(in) :: vector  ! endpoints
    TYPE(t_geographical_coordinates) , INTENT(in) :: position_local
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    REAL(wp), INTENT(out) ::  x, y
    
    CHARACTER(LEN=*), PARAMETER :: method_name='Get3DVectorToPlanarLocal'
    
    SELECT CASE(geometry_info%geometry_type)

!     CASE (planar_torus_geometry)
!       CALL finish(method_name, "planar_torus_geometry is not implemented yet")
      
    CASE (sphere_geometry)
      CALL cvec2gvec ( vector%x(1), vector%x(2), vector%x(3),     &
        & position_local%lon, position_local%lat, &
        & x, y)
    CASE ( planar_channel_geometry, planar_geometry, planar_torus_geometry )
      ! just a projection
      x = vector%x(1)
      y = vector%x(2)
      
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT
  END SUBROUTINE Get3DVectorToPlanarLocal
  !-------------------------------------------------------------------------
            
  !-------------------------------------------------------------------------
  !>
  !! returns the vector from x to y
  FUNCTION distance_vector (y, x, geometry_info)  result (d_vector)
    TYPE(t_cartesian_coordinates), INTENT(in) :: y, x  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates) :: d_vector

    REAL(wp) :: channel_x_modulo, channel_y_modulo
    CHARACTER(LEN=*), PARAMETER :: method_name='distance_vector'
    !-----------------------------------------------------------------------

    SELECT CASE(geometry_info%geometry_type)

    CASE (planar_torus_geometry)
      d_vector%x = y%x - x%x
      channel_x_modulo = geometry_info%domain_length ! / geometry_info%sphere_radius
      IF (ABS(d_vector%x(1)) > channel_x_modulo  * 0.5_wp) THEN
        d_vector%x(1) = SIGN(1.0_wp, d_vector%x(1)) * (ABS(d_vector%x(1)) - channel_x_modulo)
      ENDIF
      channel_y_modulo = geometry_info%domain_height ! / geometry_info%sphere_radius
      IF (ABS(d_vector%x(2)) > channel_y_modulo  * 0.5_wp) THEN
        d_vector%x(2) = SIGN(1.0_wp, d_vector%x(2)) * (ABS(d_vector%x(2)) - channel_y_modulo)
      ENDIF
      d_vector%x = d_vector%x / grid_sphere_radius
      
    CASE (sphere_geometry)
      d_vector%x = y%x - x%x

    CASE ( planar_channel_geometry )
      d_vector%x = y%x - x%x
      channel_x_modulo = geometry_info%domain_length ! / geometry_info%sphere_radius
      IF (ABS(d_vector%x(1)) > channel_x_modulo  * 0.5_wp) THEN
        d_vector%x(1) = SIGN(1.0_wp, d_vector%x(1)) * (ABS(d_vector%x(1)) - channel_x_modulo)
      ENDIF
      d_vector%x = d_vector%x / grid_sphere_radius
    CASE ( planar_geometry )
      d_vector%x = (y%x - x%x) / grid_sphere_radius

    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END FUNCTION distance_vector  
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! returns the distance from x to y
  REAL(wp) FUNCTION distance(y, x, geometry_info) 
    TYPE(t_cartesian_coordinates), INTENT(in) :: y, x  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates) :: d_vector

    TYPE(t_cartesian_coordinates) :: d
    !-----------------------------------------------------------------------
    d = distance_vector (y, x, geometry_info)
    distance = SQRT(SUM(d%x * d%x))

  END FUNCTION distance
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  REAL(wp) FUNCTION planar_triangle_area (x1, x2, x3, geometry_info) 
    TYPE(t_cartesian_coordinates), INTENT(in) :: x1, x2, x3  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info

    TYPE(t_cartesian_coordinates) :: d_vector_1, d_vector_2, p_vector
!     REAL(wp) :: test_planar_triangle_area
    !-----------------------------------------------------------------------
    d_vector_1 = distance_vector(x2, x1, geometry_info)
    d_vector_2 = distance_vector(x3, x1, geometry_info)
    ! the area is the one defined by the points 0, d_vector_1,  d_vector_2
    p_vector = vector_product(d_vector_1, d_vector_2)
    planar_triangle_area = 0.5_wp * SQRT(SUM(p_vector%x * p_vector%x))

    !test_planar_triangle_area  = triangle_area(x1, x2, x3)
!     write(0,*) "point 1:", x1%x
!     write(0,*) "point 2:", x2%x
!     write(0,*) "point 3:", x3%x
!     write(0,*) "vector 2-1:", d_vector_1%x
!     write(0,*) "vector 3-1:", d_vector_2%x
!     write(0,*) "triangle_area:", planar_triangle_area ! test_planar_triangle_area

  END FUNCTION planar_triangle_area
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! returns the middle between x and y
  FUNCTION planar_middle (y, x, geometry_info)  result (middle)
    TYPE(t_cartesian_coordinates), INTENT(in) :: y, x  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates) :: middle

    TYPE(t_cartesian_coordinates) :: d_vector
    !-----------------------------------------------------------------------

    d_vector = distance_vector(y, x, geometry_info)
    middle%x = x%x + d_vector%x * 0.5_wp

  END FUNCTION planar_middle
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  REAL(wp) FUNCTION planar_distance (y, x, geometry_info)  
    TYPE(t_cartesian_coordinates), INTENT(in) :: y, x  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    
    TYPE(t_cartesian_coordinates) :: d_vector

    d_vector = distance_vector(y, x, geometry_info)
    planar_distance  = SQRT(SUM(d_vector%x * d_vector%x))

  END FUNCTION planar_distance
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! returns the vector from x to y
  FUNCTION get_surface_normal (x, geometry_info)  result (z_vector)
    TYPE(t_cartesian_coordinates), INTENT(in) :: x  ! endpoints
    TYPE(t_grid_geometry_info), INTENT(in) :: geometry_info
    TYPE(t_cartesian_coordinates) :: z_vector

    CHARACTER(LEN=*), PARAMETER :: method_name='get_surface_normal'
    !-----------------------------------------------------------------------

    SELECT CASE(geometry_info%geometry_type)

    CASE (sphere_geometry)
      z_vector%x = x%x !/ SQRT(SUM(x%x * x%x))
      z_vector%x=z_vector%x/SQRT(SUM(z_vector%x * z_vector%x))
    CASE ( planar_channel_geometry, planar_geometry, planar_torus_geometry )
      z_vector%x = (/0.0_wp, 0.0_wp, 1.0_wp/)
    CASE DEFAULT
      CALL finish(method_name, "Undefined geometry type")
    END SELECT

  END FUNCTION get_surface_normal
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE construct_operators_coefficients( patch_3D, operators_coefficients, solverCoeff_sp)
    TYPE(t_patch_3D),TARGET,INTENT(inout) :: patch_3D
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    CALL allocate_operators_coefficients( patch_3d%p_patch_2d(1), operators_coefficients, solverCoeff_sp)
    CALL par_init_operator_coeff( patch_3d, operators_coefficients, solverCoeff_sp)

  END SUBROUTINE construct_operators_coefficients
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE destruct_operators_coefficients( operators_coefficients, solverCoeff_sp )
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    CALL deallocate_operators_coefficients( operators_coefficients, solverCoeff_sp )

  END SUBROUTINE destruct_operators_coefficients
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Allocation of operators coefficients.
  !
  !
!<Optimize:inUse>
  SUBROUTINE allocate_operators_coefficients( patch_2D, operators_coefficients, solverCoeff_sp)
    !
    TYPE(t_patch),TARGET,INTENT(in)       :: patch_2D
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, nz_lev
    INTEGER :: return_status,ie,i
    INTEGER :: cells_startidx, cells_endidx
    INTEGER :: edges_startidx, edges_endidx
    INTEGER :: jc,je,jk,block
    CHARACTER(len=max_char_length) :: var_suffix

    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: all_cells

    !-----------------------------------------------------------------------
    !
    ! determine size of arrays, i.e.
    ! values for the blocking
    !
    alloc_cell_blocks  = patch_2D%alloc_cell_blocks
    nblks_e  = patch_2D%nblks_e
    nblks_v  = patch_2D%nblks_v
    nz_lev   = n_zlev

    no_primal_edges = patch_2D%cells%max_connectivity
    no_dual_edges   = patch_2D%verts%max_connectivity

    ALLOCATE(operators_coefficients%div_coeff(nproma,n_zlev,alloc_cell_blocks,no_primal_edges),&
      & stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_div failed')
    ENDIF

    ALLOCATE(operators_coefficients%grad_coeff(nproma,n_zlev,nblks_e),&
      & stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for geofac_grad failed')
    ENDIF

    ALLOCATE(operators_coefficients%averageCellsToEdges(nproma,nblks_e,2),&
      & stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for averageCellsToEdges failed')
    ENDIF

    ALLOCATE(operators_coefficients%rot_coeff(nproma,n_zlev,nblks_v,no_dual_edges),&
      & stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d', &
        & 'allocation for geofac_rot failed')
    ENDIF

!     ALLOCATE(operators_coefficients%n2s_coeff(nproma,n_zlev,alloc_cell_blocks,no_primal_edges+1),&
!       & stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d',                       &
!         & 'allocation for geofac_n2s failed')
!     ENDIF
!     ALLOCATE(operators_coefficients%n2v_coeff(nproma,n_zlev,nblks_e),&
!       & stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d',                       &
!         & 'allocation for geofac_n2v failed')
!     ENDIF
    !
!     ALLOCATE(operators_coefficients%dist_cell2edge(nproma,n_zlev,nblks_e,2),stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d:allocating dist_cell2edge failed')
!     ENDIF

    ALLOCATE(operators_coefficients%vertex_bnd_edge_idx(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating vertex_bnd_edge_idx failed')
    ENDIF
    ALLOCATE(operators_coefficients%vertex_bnd_edge_blk(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating vertex_bnd_edge_blk failed')
    ENDIF
    ALLOCATE(operators_coefficients%boundaryEdge_Coefficient_Index(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge_idx failed')
    ENDIF
!     ALLOCATE(operators_coefficients%orientation(nproma,n_zlev,nblks_v,no_dual_edges-2),stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d:allocating orientation failed')
!     ENDIF
    ALLOCATE(operators_coefficients%bnd_edges_per_vertex(nproma,n_zlev,nblks_v),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating bnd_edges_per_vertex failed')
    ENDIF
!     ALLOCATE(operators_coefficients%upwind_cell_idx(nproma,n_zlev,nblks_e),stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_idx failed')
!     ENDIF
!     ALLOCATE(operators_coefficients%upwind_cell_blk(nproma,n_zlev,nblks_e),stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d:allocating upwind_cell_blk failed')
!     ENDIF
    !
    ! arrays that are required for setting up the scalar product
    !
    !coefficients for edge to cell mapping, one half of the scalar product.
    !Dimension: nproma,alloc_cell_blocks encode number of cells, 1:3 corresponds to number
    !of edges per cell, 1:2 is for u and v component of cell vector
    !     ALLOCATE(operators_coefficients%edge2cell_coeff(nproma,alloc_cell_blocks,1:3, 1:2),STAT=return_status)
    !     IF (return_status /= SUCCESS) THEN
    !       CALL finish ('allocating edge2cell_coeff failed')
    !     ENDIF
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff(nproma,nz_lev,nblks_e,1:2*no_primal_edges),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff failed')
    ENDIF

! Note: the following have the connectivity as first index
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff_top(1:2*no_primal_edges, nproma, nblks_e),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff_top failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff_integrated(1:2*no_primal_edges, nproma, nblks_e),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff_integrated failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge2edge_viacell_coeff_all(1:2*no_primal_edges, nproma, nblks_e),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2edge_viacell_coeff_all failed')
    ENDIF

    IF (select_lhs .GE. select_lhs_matrix .AND. select_lhs .LE. select_lhs_matrix + 1) THEN
      IF ( patch_2D%cells%max_connectivity /= 3 ) &
        CALL finish("select_lhs_matrix","this works only for triangles")
      ALLOCATE(operators_coefficients%lhs_all(0:9, nproma, alloc_cell_blocks),              &
               operators_coefficients%lhs_CellToCell_index(1:9, nproma, alloc_cell_blocks), &
               operators_coefficients%lhs_CellToCell_block(1:9, nproma, alloc_cell_blocks), &
                stat=return_status)
      IF (return_status /= success) &
        & CALL finish ('mo_operator_ocean_coeff_3d:allocating lhs_all failed')
    ENDIF
      
     
    
! End Note



    ALLOCATE(operators_coefficients%edge2cell_coeff_cc(nproma,nz_lev,alloc_cell_blocks,1:no_primal_edges),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

!     ALLOCATE(operators_coefficients%edge2cell_coeff_cc_dyn(nproma,1,alloc_cell_blocks,1:no_primal_edges),stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc_dyn failed')
!     ENDIF
    !ALLOCATE(operators_coefficients%edge2vert_coeff_cc_dyn(nproma,1,nblks_v,1:no_dual_edges),stat=return_status)
    !IF (return_status /= success) THEN
    !  CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff_cc_dyn failed')
    !ENDIF

    !coefficients for transposed of edge to cell mapping, second half of the scalar product.
    !Dimension: nproma,nblks_e encode number of edges, 1:2 is for cell neighbors of an edge
    ALLOCATE(operators_coefficients%edge2cell_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating transposed edge2cell_coeff failed')
    ENDIF

    !
    !coefficients for edge to vertex mapping.
    !
    !Dimension: nproma,nblks_v encode number of vertices,
    !1:6 is number of edges of a vertex,
    !1:2 is for u and v component of vertex vector
    ALLOCATE(operators_coefficients%edge2vert_coeff_cc(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF

    ALLOCATE(operators_coefficients%edge2vert_coeff_cc_t(nproma,nz_lev,nblks_e,1:2),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_coeff failed')
    ENDIF
    ALLOCATE(operators_coefficients%edge2vert_vector_cc(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2vert_vector failed')
    ENDIF

   ALLOCATE(operators_coefficients%edge2edge_viavert_coeff(nproma,nz_lev,nblks_e,1:2*no_dual_edges),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating edge2cell_coeff_cc failed')
    ENDIF

    !
    !normalizing factors for edge to cell mapping.
    !
    !Either by fixed volume or by variable one taking the surface elevation
    !into account. The later one depends on time and space.
    ALLOCATE(operators_coefficients%fixed_vol_norm(nproma,nz_lev,alloc_cell_blocks),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating fixed_vol_norm failed')
    ENDIF
!     ALLOCATE(operators_coefficients%variable_vol_norm(nproma,nz_lev,alloc_cell_blocks,1:no_primal_edges),stat=return_status)
!     IF (return_status /= success) THEN
!     ENDIF

!     ALLOCATE(operators_coefficients%variable_dual_vol_norm(nproma,nz_lev,nblks_v,1:no_dual_edges),stat=return_status)
!     IF (return_status /= success) THEN
!       CALL finish ('mo_operator_ocean_coeff_3d:allocating variable_dual_vol_norm failed')
!     ENDIF

    !---------------------------------------------------------------
    ! allocate t_verticalAdvection_ppm_coefficients
    ALLOCATE(operators_coefficients%verticalAdvectionPPMcoeffs(alloc_cell_blocks),stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d:allocating verticalAdvectionPPMcoeffs failed')
    ENDIF


!ICON_OMP_PARALLEL_DO PRIVATE(return_status) ICON_OMP_DEFAULT_SCHEDULE
    DO block = 1, alloc_cell_blocks
      ALLOCATE( &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_This_toBelow(nproma,nz_lev),                &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_This_toThisBelow(nproma,nz_lev),            &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeight_2xBelow_x_RatioThis_toThisBelow(nproma,nz_lev),  &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_This_toThisAboveBelow(nproma,nz_lev),       &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_2xAboveplusThis_toThisBelow(nproma,nz_lev), &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_2xBelowplusThis_toThisAbove(nproma,nz_lev), &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_ThisAbove_to2xThisplusBelow(nproma,nz_lev), &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_ThisBelow_to2xThisplusAbove(nproma,nz_lev), &
        & operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeight_inv_ThisAboveBelow2Below(nproma,nz_lev),         &
        & stat=return_status)
      IF (return_status /= success) THEN
        CALL finish ('mo_operator_ocean_coeff_3d:allocating verticalAdvectionPPMcoeffs failed')
      ENDIF
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_This_toBelow(:,:)                = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_This_toThisBelow(:,:)            = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeight_2xBelow_x_RatioThis_toThisBelow(:,:)  = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_This_toThisAboveBelow(:,:)       = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_2xAboveplusThis_toThisBelow(:,:) = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_2xBelowplusThis_toThisAbove(:,:) = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_ThisAbove_to2xThisplusBelow(:,:) = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeightRatio_ThisBelow_to2xThisplusAbove(:,:) = 0.0_wp
      operators_coefficients%verticalAdvectionPPMcoeffs(block)%cellHeight_inv_ThisAboveBelow2Below(:,:)         = 0.0_wp
      
    ENDDO
!ICON_OMP_END_PARALLEL_DO
    


    !---------------------------------------------------------------
    ! allocate single precision operators
    IF (select_solver == select_gmres_mp_r) THEN

      ALLOCATE(solverCoeff_sp%div_coeff(nproma,alloc_cell_blocks,no_primal_edges),&
        & solverCoeff_sp%grad_coeff(nproma,nblks_e),                    &
        & solverCoeff_sp%edge2edge_viacell_coeff_all(1:2*no_primal_edges, nproma, nblks_e), &
        & solverCoeff_sp%edge_thickness(nproma, nblks_e),                      &
        & solverCoeff_sp%cell_thickness(nproma, alloc_cell_blocks),                      &
        & stat=return_status)

      IF (return_status /= success) THEN
        CALL finish ('mo_operator_ocean_coeff_3d',                 &
          & 'allocation for solverCoeff_sp failed')
      ENDIF


    ENDIF

    ALLOCATE(operators_coefficients%edges_SeaBoundaryLevel(nproma,n_zlev,nblks_e),&
      & operators_coefficients%cells_SeaBoundaryLevel(nproma,n_zlev,alloc_cell_blocks), stat=return_status)
    IF (return_status /= success) THEN
      CALL finish ('mo_operator_ocean_coeff_3d',                 &
        & 'allocation for edges_SeaBoundaryLevel failed')
    ENDIF

    !---------------------------------------------------------------
    !
    ! initialize all components
    !
    DO ie = 1,3
      operators_coefficients%edge2cell_coeff_cc%x(ie)     = 0._wp
      operators_coefficients%edge2cell_coeff_cc_t%x(ie)   = 0._wp
      operators_coefficients%edge2vert_coeff_cc%x(ie)     = 0._wp
      operators_coefficients%edge2vert_coeff_cc_t%x(ie)   = 0._wp
      operators_coefficients%edge2vert_vector_cc%x(ie)    = 0._wp
!       operators_coefficients%edge2cell_coeff_cc_dyn%x(ie) = 0._wp
      !operators_coefficients%edge2vert_coeff_cc_dyn%x(ie) = 0._wp
    END DO

    all_cells => patch_2D%cells%all
    all_edges => patch_2D%edges%all
    !all_verts => patch_2D%verts%all

    operators_coefficients%edge2edge_viacell_coeff= 0._wp
    operators_coefficients%edge2edge_viavert_coeff= 0._wp

    operators_coefficients%fixed_vol_norm         = 0._wp
!     operators_coefficients%variable_vol_norm      = 0._wp
!     operators_coefficients%variable_dual_vol_norm = 0._wp

!     operators_coefficients%dist_cell2edge = 0._wp

    operators_coefficients%div_coeff  = 0._wp
    operators_coefficients%rot_coeff  = 0._wp
    operators_coefficients%grad_coeff = 0._wp
    operators_coefficients%averageCellsToEdges = 0._wp

    operators_coefficients%vertex_bnd_edge_idx = 0
    operators_coefficients%vertex_bnd_edge_blk = 0
    operators_coefficients%boundaryEdge_Coefficient_Index     = 0
!     operators_coefficients%orientation  = 0.0_wp
    operators_coefficients%bnd_edges_per_vertex= 0

!    operators_coefficients%upwind_cell_idx = 1
!    operators_coefficients%upwind_cell_blk = 1

    CALL message ('mo_operator_ocean_coeff_3d:allocate_operators_coefficients',&
      & 'memory allocation finished')

  END SUBROUTINE allocate_operators_coefficients
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Deallocation of operators coefficients.
!<Optimize:inUse>
  SUBROUTINE deallocate_operators_coefficients( operators_coefficients, solverCoeff_sp )
    ! !
    TYPE(t_operator_coeff), INTENT(inout) :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp

    !-----------------------------------------------------------------------
    DEALLOCATE(operators_coefficients%div_coeff)
    DEALLOCATE(operators_coefficients%grad_coeff)
    DEALLOCATE(operators_coefficients%averageCellsToEdges)

    DEALLOCATE(operators_coefficients%rot_coeff)

!     DEALLOCATE(operators_coefficients%n2s_coeff)
!     DEALLOCATE(operators_coefficients%n2v_coeff)
    !
!     DEALLOCATE(operators_coefficients%dist_cell2edge)

    DEALLOCATE(operators_coefficients%vertex_bnd_edge_idx)
    DEALLOCATE(operators_coefficients%vertex_bnd_edge_blk)
    DEALLOCATE(operators_coefficients%boundaryEdge_Coefficient_Index)
!     DEALLOCATE(operators_coefficients%orientation)
    DEALLOCATE(operators_coefficients%bnd_edges_per_vertex)
!     DEALLOCATE(operators_coefficients%upwind_cell_idx)
!     DEALLOCATE(operators_coefficients%upwind_cell_blk)

    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff)
    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff_top)
    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff_integrated)
    DEALLOCATE(operators_coefficients%edge2edge_viacell_coeff_all)
    IF (select_lhs .GE. select_lhs_matrix .AND. select_lhs .LE. select_lhs_matrix + 1) &
      & DEALLOCATE(operators_coefficients%lhs_all, &
        & operators_coefficients%lhs_CellToCell_index, &
        & operators_coefficients%lhs_CellToCell_block)
 
    DEALLOCATE(operators_coefficients%edge2cell_coeff_cc)

!     DEALLOCATE(operators_coefficients%edge2cell_coeff_cc_dyn)

    DEALLOCATE(operators_coefficients%edge2cell_coeff_cc_t)

    DEALLOCATE(operators_coefficients%edge2vert_coeff_cc)

    DEALLOCATE(operators_coefficients%edge2vert_coeff_cc_t)
    DEALLOCATE(operators_coefficients%edge2vert_vector_cc)

    DEALLOCATE(operators_coefficients%edge2edge_viavert_coeff)

    DEALLOCATE(operators_coefficients%fixed_vol_norm)
!     DEALLOCATE(operators_coefficients%variable_vol_norm)
!     DEALLOCATE(operators_coefficients%variable_dual_vol_norm)

    IF (select_solver == select_gmres_mp_r) THEN

      DEALLOCATE(solverCoeff_sp%div_coeff,              &
        & solverCoeff_sp%grad_coeff,                    &
        & solverCoeff_sp%edge2edge_viacell_coeff_all,   &
        & solverCoeff_sp%edge_thickness,                &
        & solverCoeff_sp%cell_thickness)

    ENDIF

    DEALLOCATE(operators_coefficients%edges_SeaBoundaryLevel)
    DEALLOCATE(operators_coefficients%cells_SeaBoundaryLevel)

  END SUBROUTINE deallocate_operators_coefficients
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !> Initialize expansion coefficients.
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE par_init_operator_coeff( patch_3D, operators_coefficients, solverCoeff_sp)
    !
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: patch_3D
    TYPE(t_operator_coeff),   INTENT(inout) :: operators_coefficients
    TYPE(t_solverCoeff_singlePrecision), INTENT(inout) :: solverCoeff_sp
   !
   !Local variables:
   TYPE(t_patch),POINTER :: patch_2D
    !-----------------------------------------------------------------------
    !TYPE(t_cartesian_coordinates) :: check_v(nproma, n_zlev, patch_2D%nblks_v, 6)
    !REAL(wp) :: check_r(nproma, n_zlev, patch_2D%alloc_cell_blocks, 3)
    !REAL(wp) :: max_diff, max_val
    !-----------------------------------------------------------------------
    patch_2D => patch_3D%p_patch_2D(1)

    CALL init_geometry_forCoefficients(patch_2D)

    CALL init_operator_coeffs( patch_2D, operators_coefficients)

    CALL init_diff_operator_coeff_3D ( patch_2D, operators_coefficients )

    CALL apply_boundary2coeffs(patch_3D, operators_coefficients)

    CALL init_verticalAdvection_ppm_coefficients(patch_3D, operators_coefficients%verticalAdvectionPPMcoeffs)
    
    IF (select_solver == select_gmres_mp_r) THEN
      solverCoeff_sp%grad_coeff(:,:)   = REAL(operators_coefficients%grad_coeff(:,1,:), sp)
      solverCoeff_sp%div_coeff(:,:,:)  = REAL(operators_coefficients%div_coeff(:,1,:,:), sp)
      solverCoeff_sp%edge2edge_viacell_coeff_all(:,:,:)  = REAL(operators_coefficients%edge2edge_viacell_coeff_all(:,:,:), sp)
    ENDIF

    CALL destruct_geometry_forCoefficients()

  END SUBROUTINE par_init_operator_coeff
  !-------------------------------------------------------------------------

  !---------------------------------------------------------------------------------
  !>
  !!
!<Optimize:inUse>
  SUBROUTINE init_verticalAdvection_ppm_coefficients( patch_3D, vertAdvPPM_coeffs)
    !
    ! patch_2D on which computation is performed
    TYPE(t_patch_3D ),TARGET, INTENT(IN)   :: patch_3D
    TYPE(t_verticalAdvection_ppm_coefficients), TARGET :: vertAdvPPM_coeffs(:)

    !  local variables
    INTEGER            :: cell_StartIndex, cell_EndIndex
    INTEGER            :: jc, jb, je, jk
    INTEGER            :: thisLevel, levelAbove, levelBelow, level2Below, cell_levels

    TYPE(t_patch), POINTER  :: patch_2D
    TYPE(t_subset_range), POINTER :: all_cells

    REAL(wp), POINTER :: cell_thickeness(:,:,:)
    !-------------------------------------------------------------------------------
    ! pointers for the ppm vertical transport
    TYPE(t_verticalAdvection_ppm_coefficients), POINTER :: vertAdvPPM
    !-------------------------------------------------------------------------------
    !CALL message (TRIM(routine), 'start')
    patch_2D            => patch_3D%p_patch_2D(1)
    all_cells           => patch_2D%cells%all
    cell_thickeness     => patch_3D%p_patch_1d(1)%prism_thick_c

    !-------------------------------------------------------------------------
    ! update the coefficients for the upwind_vflux_ppm_fast vertical advection
!ICON_OMP_PARALLEL_DO PRIVATE(cell_StartIndex, cell_EndIndex, vertAdvPPM, jc, cell_levels, thisLevel, levelAbove, &
!ICON_OMP  levelBelow, level2Below   ) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, cell_StartIndex, cell_EndIndex)
      vertAdvPPM => vertAdvPPM_coeffs(jb)
      DO jc = cell_StartIndex, cell_EndIndex

        cell_levels = patch_3D%p_patch_1D(1)%dolic_c(jc,jb)

        DO thisLevel  = 1, cell_levels-1
          levelBelow = thisLevel + 1

          vertAdvPPM%cellHeightRatio_This_toBelow(jc, thisLevel) = &
            & cell_thickeness(jc, thisLevel, jb) / cell_thickeness(jc, levelBelow, jb)

          vertAdvPPM%cellHeightRatio_This_toThisBelow(jc, thisLevel) = &
            & cell_thickeness(jc, thisLevel, jb) / &
            & (cell_thickeness(jc, thisLevel, jb) + cell_thickeness(jc, levelBelow, jb))

          vertAdvPPM%cellHeight_2xBelow_x_RatioThis_toThisBelow(jc,thisLevel) = &
            & 2._wp * cell_thickeness(jc,levelBelow, jb) * &
            & vertAdvPPM%cellHeightRatio_This_toThisBelow(jc, thisLevel)

        ENDDO


        DO thisLevel  = 2, cell_levels-1
          levelAbove = thisLevel - 1
          levelBelow = thisLevel + 1


          vertAdvPPM%cellHeightRatio_This_toThisAboveBelow(jc,thisLevel) = &
            & cell_thickeness(jc, thisLevel ,jb) / &
            &   (cell_thickeness(jc,levelAbove,jb) + cell_thickeness(jc,thisLevel,jb)    &
            &    + cell_thickeness(jc,levelBelow,jb))

          vertAdvPPM%cellHeightRatio_2xAboveplusThis_toThisBelow(jc,thisLevel) = &
            & (2._wp * cell_thickeness(jc,levelAbove,jb) + cell_thickeness(jc,thisLevel,jb))     &
            & / (cell_thickeness(jc,levelBelow,jb) + cell_thickeness(jc,thisLevel,jb))

          vertAdvPPM%cellHeightRatio_2xBelowplusThis_toThisAbove(jc,thisLevel) = &
            & + (cell_thickeness(jc,thisLevel,jb) + 2._wp * cell_thickeness(jc,levelBelow,jb))   &
            & / (cell_thickeness(jc,levelAbove,jb) + cell_thickeness(jc,thisLevel,jb))

          vertAdvPPM%cellHeightRatio_ThisAbove_to2xThisplusBelow(jc,thisLevel) =                         &
            &  (cell_thickeness(jc,levelAbove,jb) + cell_thickeness(jc,thisLevel,jb))            &
            & / (2._wp*cell_thickeness(jc,thisLevel,jb) + cell_thickeness(jc,levelBelow,jb))

          vertAdvPPM%cellHeightRatio_ThisBelow_to2xThisplusAbove(jc,thisLevel) =                 &
            &  (cell_thickeness(jc,levelBelow,jb) + cell_thickeness(jc,thisLevel,jb))                  &
            & / (2._wp*cell_thickeness(jc,thisLevel,jb) + cell_thickeness(jc,levelAbove,jb))
            ! = 1 / cellHeightRatio_2xBelowplusThis_toThisAbove(levelBelow)

        ENDDO

        DO thisLevel  = 2, cell_levels-2
          levelAbove  = thisLevel - 1
          levelBelow  = thisLevel + 1
          level2Below = thisLevel + 2

          vertAdvPPM%cellHeight_inv_ThisAboveBelow2Below(jc,thisLevel) =                                  &
            & 1._wp / (cell_thickeness(jc,levelAbove,jb) + cell_thickeness(jc,thisLevel,jb)       &
            &           + cell_thickeness(jc,levelBelow,jb) + cell_thickeness(jc,level2Below,jb))

        ENDDO

      END DO
    END DO
!ICON_OMP_END_PARALLEL_DO

  END SUBROUTINE init_verticalAdvection_ppm_coefficients
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE init_operator_coeffs( patch_2D, operators_coefficients)
    TYPE(t_patch)    , TARGET, INTENT(INOUT)     :: patch_2D
    TYPE(t_operator_coeff),    INTENT(inout)     :: operators_coefficients

!Local variables

    REAL(wp)                      :: div_coeff              (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)
    REAL(wp)                      :: rot_coeff              (1:nproma,1:patch_2D%nblks_v,1:no_dual_edges)
    REAL(wp)                      :: grad_coeff             (1:nproma,1:patch_2D%nblks_e)

    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc     (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)


    TYPE(t_subset_range), POINTER :: owned_edges         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_cells         ! these are the owned entities
    TYPE(t_subset_range), POINTER :: owned_verts         ! these are the owned entities
    TYPE(t_cartesian_coordinates) ::  edge_center !, vertex_center
    TYPE(t_cartesian_coordinates) :: dist_vector!, dist_vector_basic
    TYPE(t_cartesian_coordinates) :: coriolis_cartesian_coordinates
    TYPE(t_geographical_coordinates) :: coriolis_geo_coordinates, geo_coordinates
    REAL(wp) :: basin_center_lat_rad, basin_height_rad
    REAL(wp) :: length
    REAL(wp) :: inverse_sphere_radius
    REAL(wp) :: w1, w2

    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: cell_1_index, cell_1_block, cell_2_index, cell_2_block
    INTEGER :: vertex_1_index, vertex_1_block, vertex_2_index, vertex_2_block
    INTEGER :: level
    CHARACTER(*), PARAMETER :: method_name = "init_operator_coeffs"
    !-----------------------------------------------------------------------
    inverse_sphere_radius = 1.0_wp / grid_sphere_radius

    owned_edges => patch_2D%edges%owned
    owned_cells => patch_2D%cells%owned
    owned_verts => patch_2D%verts%owned

    rot_coeff(:,:,:)        = 0.0_wp
    div_coeff(:,:,:)        = 0.0_wp
    grad_coeff(:,:)         = 0.0_wp
    !-------------------------------------------



    !-------------------------------------------
    !2) calculate coefficients for difference operators
    !
    !2a) divergence
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        DO neigbor=1, patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges

          edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x = 0.0_wp

          edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

          div_coeff(cell_index,cell_block,neigbor) =                           &
              & patch_2D%edges%primal_edge_length(edge_index,edge_block) *        &
              & patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)  / &
              & patch_2D%cells%area(cell_index,cell_block)

        ENDDO !neigbor=1,patch_2D%num_edges
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block
 

   !2b) gradient, average
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        grad_coeff(edge_index,edge_block)&
          & =1.0_wp/ dual_edge_length(edge_index,edge_block)!patch_2D%edges%inv_dual_edge_length(edge_index, edge_block)

        w1 = 0.0_wp
        w2 = 0.0_wp
        IF (dist_cell2edge(edge_index,edge_block,1) > 0.0_wp) THEN
          w1 = exp(-dist_cell2edge(edge_index,edge_block,1))
        ENDIF
        IF (dist_cell2edge(edge_index,edge_block,2) > 0.0_wp) THEN
          w2 = exp(-dist_cell2edge(edge_index,edge_block,2))
        ENDIF
         
        operators_coefficients%averageCellsToEdges(edge_index,edge_block,1) = w1 / (w1+w2)
        operators_coefficients%averageCellsToEdges(edge_index,edge_block,2) = w2 / (w1+w2)

      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    
!    CALL dbg_print('dual_edge_length',dual_edge_length, &
!      & this_mod_name,1, in_subset=owned_edges)
   !2c) curl coefficients
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        DO neigbor=1, patch_2D%verts%num_edges(vertex_index,vertex_block)!no_dual_edges

          edge_index = patch_2D%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch_2D%verts%edge_blk(vertex_index, vertex_block, neigbor)

          IF (edge_block > 0) THEN
            rot_coeff(vertex_index,vertex_block,neigbor)           &
             = dual_edge_length(edge_index,edge_block) * grid_sphere_radius &     !PK why not dual_edge_length from above ??
!             & = patch_2D%edges%primal_edge_length(edge_index,edge_block) &         !PK why not dual_edge_length from above ??
             &  * patch_2D%verts%edge_orientation(vertex_index,vertex_block,neigbor)

            dist_vector = distance_vector( &
              & patch_2D%edges%cartesian_center(edge_index,edge_block), &
              & patch_2D%verts%cartesian(vertex_index, vertex_block), &
              & patch_2D%geometry_info)
            IF (DOT_PRODUCT(patch_2D%edges%dual_cart_normal(edge_index,edge_block)%x, dist_vector%x) &
                 * patch_2D%verts%edge_orientation(vertex_index,vertex_block,neigbor) < 0.0_wp ) THEN
              write(0,*) "dist_vector: ", dist_vector%x
              write(0,*) "dual_cart_normal: ", patch_2D%edges%dual_cart_normal(edge_index,edge_block)%x
              write(0,*) "orientation: ",  patch_2D%verts%edge_orientation(vertex_index,vertex_block,neigbor)
              write(0,*) "DOT_PRODUCT: ", DOT_PRODUCT(patch_2D%edges%dual_cart_normal(edge_index,edge_block)%x, dist_vector%x)
              CALL finish(method_name, "wrong orientation for rot_coeff")
            ENDIF

          ENDIF
        ENDDO !neigbor=1,6
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
   
    !Copy coefficients to 3D
    DO level=1,n_zlev
      operators_coefficients%div_coeff(:,level,:,:) = div_coeff(:,:,:)
      operators_coefficients%rot_coeff(:,level,:,:) = rot_coeff(:,:,:)
      operators_coefficients%grad_coeff(:,level,:)  = grad_coeff(:,:)
    END DO
    !-------------------
    ! sync the results
    CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%grad_coeff(:,:,:))
    CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%averageCellsToEdges(:,:,1))
    CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%averageCellsToEdges(:,:,2))
    DO neigbor=1,no_primal_edges
      CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%div_coeff(:,:,:,neigbor))
    END DO
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%rot_coeff(:,:,:,neigbor))
    END DO

    !----------------------------------------------------
    CALL init_operator_coeffs_cell( patch_2D, operators_coefficients )
    CALL init_operator_coeffs_vertex( patch_2D, operators_coefficients)

    !----------------------------------------------------
    ! 9) recalculate the coriolis coefficient
    ! It is required if we use the middle of the dual_edge_length
!     IF (MID_POINT_DUAL_EDGE) THEN
! 
!       IF (CORIOLIS_TYPE == full_coriolis) THEN
! 
!         DO edge_block = owned_edges%start_block, owned_edges%end_block
!           CALL get_index_range(owned_edges, edge_block, start_index, end_index)
!           DO edge_index = start_index, end_index
! 
!              coriolis_geo_coordinates = cc2gc(dual_edge_middle(edge_index,edge_block))
!              patch_2D%edges%f_e(edge_index,edge_block) = &
!                & 2._wp * grid_angular_velocity * SIN(coriolis_geo_coordinates%lat)
! 
!           ENDDO
!         ENDDO
! 
!       ELSEIF (CORIOLIS_TYPE == BETA_PLANE_CORIOLIS) THEN
! 
!         basin_center_lat_rad = basin_center_lat * deg2rad
!         basin_height_rad     = basin_height_deg * deg2rad
!         coriolis_geo_coordinates%lat = basin_center_lat_rad - 0.5_wp * basin_height_rad
!         coriolis_geo_coordinates%lon = 0.0_wp
!         coriolis_cartesian_coordinates  = gc2cc(coriolis_geo_coordinates)
! 
!         DO edge_block = owned_edges%start_block, owned_edges%end_block
!           CALL get_index_range(owned_edges, edge_block, start_index, end_index)
!           DO edge_index = start_index, end_index
! 
!           geo_coordinates     = cc2gc(dual_edge_middle(edge_index,edge_block))
!           geo_coordinates%lon = 0.0_wp
!           edge_center         = gc2cc(geo_coordinates)
!           length              = grid_sphere_radius * &
!             & arc_length(edge_center, coriolis_cartesian_coordinates)
! 
!           patch_2D%edges%f_e(edge_index,edge_block) =  2.0_wp * grid_angular_velocity * &
!             & ( sin(basin_center_lat_rad) + (cos(basin_center_lat_rad) / &
!             &   grid_sphere_radius) * length)
! 
!           ENDDO
!         ENDDO
! 
!       ENDIF !(CORIOLIS_TYPE==full_coriolis)
!     ENDIF ! (MID_POINT_DUAL_EDGE)
!     !-------------------
!     ! sync patch_2D%edges%f_e
!     CALL sync_patch_array(SYNC_E, patch_2D, patch_2D%edges%f_e)
!     !---------------------------------------------------------

  END SUBROUTINE init_operator_coeffs
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
!<Optimize:inUse>
  SUBROUTINE init_operator_coeffs_cell( patch_2D, operators_coefficients )
    TYPE(t_patch), TARGET,  INTENT(INOUT)  :: patch_2D
    TYPE(t_operator_coeff), INTENT(inout)  :: operators_coefficients

!Local variables
!
    REAL(wp)                      :: edge2edge_viacell_coeff_2D(1:nproma,1:patch_2D%nblks_e,1:2*no_primal_edges)
    !REAL(wp)                      :: dist_cell2edge         (1:nproma,1:patch_2D%nblks_e,1:2)
    REAL(wp)                      :: fixed_vol_norm         (1:nproma,patch_2D%alloc_cell_blocks)
!     REAL(wp)                      :: variable_vol_norm      (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)
    REAL(wp)                      :: norm, orientation
    REAL(wp)                      :: dist_edge_cell, dist_edge_cell_basic

    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc     (1:nproma,1:patch_2D%alloc_cell_blocks,1:no_primal_edges)
    TYPE(t_cartesian_coordinates) :: edge2cell_coeff_cc_t   (1:nproma,1:patch_2D%nblks_e,1:2)
    TYPE(t_cartesian_coordinates) :: cell_center, edge_center
    TYPE(t_cartesian_coordinates) :: dist_vector, dist_vector_basic

    INTEGER :: edge_block_cell, edge_index_cell, ictr
    INTEGER :: cell_edge
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: level

    TYPE(t_subset_range), POINTER :: owned_edges, all_edges         
    TYPE(t_subset_range), POINTER :: owned_cells, all_cells        
    CHARACTER(*), PARAMETER :: method_name = "init_operator_coeffs_cell"
    !-----------------------------------------------------------------------
    owned_edges => patch_2D%edges%owned
    all_edges   => patch_2D%edges%all
    owned_cells => patch_2D%cells%owned
    all_cells => patch_2D%cells%all
    !-------------------------------------------
    ! 3) compute:
    !   edge2cell_coeff_cc
    !   fixed_vol_norm
    !   variable_vol_norm
    edge2cell_coeff_cc(:,:,:)%x(1) = 0.0_wp
    edge2cell_coeff_cc(:,:,:)%x(2) = 0.0_wp
    edge2cell_coeff_cc(:,:,:)%x(3) = 0.0_wp

    edge2cell_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
    edge2cell_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
    edge2cell_coeff_cc_t(:,:,:)%x(3) = 0.0_wp

    fixed_vol_norm(:,:)       = 0.0_wp
!     variable_vol_norm(:,:,:)  = 0.0_wp
    edge2edge_viacell_coeff_2D(:,:,:) = 0.0_wp

    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index

        cell_center%x = patch_2D%cells%cartesian_center(cell_index, cell_block)%x
        fixed_vol_norm(cell_index,cell_block) = 0.0_wp

        !-------------------------------
        DO neigbor=1, patch_2D%cells%num_edges(cell_index,cell_block)!no_primal_edges

          edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x = 0.0_wp
!           variable_vol_norm(cell_index, cell_block, neigbor) =  0.0_wp

          edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

          IF (edge_block > 0 ) THEN
            ! we have an edge
            dist_vector = distance_vector( &
              & patch_2D%edges%cartesian_center(edge_index,edge_block), &
              & cell_center, &
              & patch_2D%geometry_info)

            norm  = SQRT(SUM( dist_vector%x * dist_vector%x))
            ! compute edge2cell_coeff_cc
            edge2cell_coeff_cc(cell_index,cell_block,neigbor)%x =  &
              & dist_vector%x *                                             &
              & prime_edge_length(edge_index,edge_block) *                  &
              & patch_2D%cells%edge_orientation(cell_index,cell_block,neigbor)

            fixed_vol_norm(cell_index,cell_block) = &
              & fixed_vol_norm(cell_index,cell_block) + &
              & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)



!             variable_vol_norm(cell_index, cell_block, neigbor) = &
!               & 0.5_wp * norm * prime_edge_length(edge_index,edge_block)

          ENDIF !(edge_block > 0 )
        ENDDO !neigbor=1,patch_2D%num_edges
        !-------------------------------
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block = owned_cells%start_block, owned_cells%end_block


    !-------------------
    ! sync the results
    CALL sync_patch_array(SYNC_C, patch_2D, fixed_vol_norm(:,:))
    DO neigbor=1,patch_2D%cells%max_connectivity
      CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_C, patch_2D, edge2cell_coeff_cc(:,:,neigbor)%x(3))
!       CALL sync_patch_array(SYNC_C, patch_2D, variable_vol_norm(:,:,neigbor))
    ENDDO
    !-------------------

   !copy calculated 2D arrays to 3D structure
    DO cell_block = all_cells%start_block, all_cells%end_block
      DO level = 1, n_zlev

       operators_coefficients%fixed_vol_norm(:,level,cell_block) = fixed_vol_norm(:,cell_block)

       DO neigbor=1,patch_2D%cells%max_connectivity

         operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(1)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(1)

         operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(2)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(2)

         operators_coefficients%edge2cell_coeff_cc(:,level,cell_block,neigbor)%x(3)  &
           &= edge2cell_coeff_cc(:,cell_block,neigbor)%x(3)

!          operators_coefficients%variable_vol_norm(:,level,cell_block,neigbor)  &
!            &= variable_vol_norm(:,cell_block,neigbor)

        ENDDO ! neigbor=1,patch_2D%cells%max_connectivity
      ENDDO  !  level = 1, n_zlev
    ENDDO ! cell_block
! no need for sync
!     CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%fixed_vol_norm(:,:,:))
!     DO neigbor=1,no_primal_edges
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%edge2cell_coeff_cc(:,:,:,neigbor)%x(1))
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%edge2cell_coeff_cc(:,:,:,neigbor)%x(2))
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%edge2cell_coeff_cc(:,:,:,neigbor)%x(3))
!       CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%variable_vol_norm(:,:,:,neigbor))
!     ENDDO
   ! output print level (1-5, fix)
   idt_src=5  
   CALL dbg_print('scalarprod: fixed_vol_norm',operators_coefficients%fixed_vol_norm, &
     & this_mod_name,idt_src, in_subset=all_cells)
    !-------------------------------------------
    ! 4) compute:
    !   edge2cell_coeff_cc_t

    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = patch_2D%edges%cartesian_center(edge_index, edge_block)%x

        DO neigbor=1,2

          cell_index = patch_2D%edges%cell_idx(edge_index, edge_block, neigbor)
          cell_block = patch_2D%edges%cell_blk(edge_index, edge_block, neigbor)

          IF (cell_index > 0) THEN

            dist_vector =  distance_vector( edge_center,                &
              & patch_2D%cells%cartesian_center(cell_index, cell_block),  &
              & patch_2D%geometry_info)

            orientation = DOT_PRODUCT(dist_vector%x, &
              & patch_2D%edges%primal_cart_normal(edge_index, edge_block)%x)              
            IF (orientation < 0.0_wp) dist_vector%x = - dist_vector%x
            IF (orientation * (1.5_wp - REAL(neigbor, wp)) <=0) &
              & CALL finish(method_name, "wrong orientation in edge2cell_coeff_cc_t")

!             edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x = &
!               dist_vector%x / dual_edge_length(edge_index, edge_block)

            edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x = &
              & patch_2D%edges%primal_cart_normal(edge_index, edge_block)%x &
              & * SQRT(SUM(dist_vector%x * dist_vector%x)) &
              & / dual_edge_length(edge_index, edge_block)

          ENDIF ! (cell_block > 0)
        ENDDO ! neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2cell_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    !   edge2cell_coeff_cc_t is computed

    !copy 2D to 3D structure
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index =  start_index, end_index
        DO level = 1, n_zlev

          operators_coefficients%edge2cell_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,1)%x

          operators_coefficients%edge2cell_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &= edge2cell_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2cell_coeff_cc_t(:,:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    !-------------------------------------------


    !-------------------------------------------
    ! 5) compute
    !   calculate edge2edge_viacell_coeff_2D
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = patch_2D%edges%cartesian_center(edge_index, edge_block)%x

        !ictr=0
        DO neigbor=1,2

          cell_index    = patch_2D%edges%cell_idx(edge_index, edge_block, neigbor)
          cell_block    = patch_2D%edges%cell_blk(edge_index, edge_block, neigbor)

          IF (cell_index <= 0) CYCLE

          IF(neigbor==1) ictr = 0
          IF(neigbor==2) ictr = no_primal_edges

          cell_center%x = patch_2D%cells%cartesian_center(cell_index, cell_block)%x

          !dist_vector_basic%x = edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x
          dist_vector_basic = distance_vector( edge_center, cell_center, patch_2D%geometry_info)

          dist_edge_cell_basic  = SQRT(SUM( dist_vector_basic%x * dist_vector_basic%x))
          dist_vector_basic%x = dist_vector_basic%x / dist_edge_cell_basic

          orientation = DOT_PRODUCT(dist_vector_basic%x, &
          & patch_2D%edges%primal_cart_normal(edge_index,edge_block)%x)
          IF (orientation < 0.0_wp) dist_vector_basic%x = - dist_vector_basic%x
          IF (orientation * (1.5_wp - REAL(neigbor, wp)) <=0) &
            & CALL finish(method_name, "wrong orientation in edge2edge_viacell_coeff_2D")
          !loop over the edges of neighbor 1 and 2
          DO cell_edge=1,patch_2D%cells%num_edges(cell_index,cell_block)

            ictr=ictr+1
            !actual edge
            edge_index_cell = patch_2D%cells%edge_idx(cell_index, cell_block, cell_edge)
            edge_block_cell = patch_2D%cells%edge_blk(cell_index, cell_block, cell_edge)

            !dist_vector%x = edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x
            dist_vector =  distance_vector( &
              & patch_2D%edges%cartesian_center(edge_index_cell, edge_block_cell),  &
              & cell_center,                &
              & patch_2D%geometry_info)

            dist_edge_cell  = SQRT(SUM( dist_vector%x * dist_vector%x))
            dist_vector%x = dist_vector%x / dist_edge_cell
            dist_vector%x = dist_vector%x * patch_2D%cells%edge_orientation(cell_index,cell_block,cell_edge)

            !This is the cosine of the angle between vectors from cell center
            !to cell edges
            edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)&
            & =DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)

            !IF(abs(edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)-1.0_wp)<1.0E-6_wp)THEN
            !  write(*,*)'ran into'
              !edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)=1.0_wp
            !ENDIF

            !multiply the cosine by length and orientation and divide by
            !dual length
            edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)=       &
              &  edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr)   &
              &* prime_edge_length(edge_index_cell,edge_block_cell)       &
              &* dist_edge_cell * dist_edge_cell_basic                    &
              &/ dual_edge_length(edge_index, edge_block)

! IF(edge_index==1.and.edge_block==1)THEN
! write(123,*)'actual angle',neigbor, edge_index_cell, edge_block_cell,ictr,&
! & edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr),&
! &DOT_PRODUCT(edge2cell_coeff_cc(cell_index,cell_block,cell_edge)%x,edge2cell_coeff_cc_t(edge_index, edge_block, neigbor)%x),&
! &acos(edge2edge_viacell_coeff_2D(edge_index,edge_block,ictr))*rad2deg
! !IF(edge_index_cell==edge_index.and.edge_block_cell==edge_block)THEN
! !write(123,*)'vecs',neigbor,dist_vector_basic%x,dist_vector%x
! !ENDIF
! ENDIF
          END DO
        ENDDO ! neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

    ! these coeffecients will not be used for-non owened edges,
    ! sync only for safety
    DO ictr=1, 2*no_primal_edges
      CALL sync_patch_array(SYNC_E, patch_2D, edge2edge_viacell_coeff_2D(:,:,ictr))
    ENDDO

    !copy 2D to 3D structure
    DO edge_block = all_edges%start_block, all_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(all_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          operators_coefficients%edge2edge_viacell_coeff(edge_index,level,edge_block,1:2*no_primal_edges) = &
            & edge2edge_viacell_coeff_2D(edge_index,edge_block,1:2*no_primal_edges)

        ENDDO
      ENDDO
    ENDDO    
  !-------------------------------------------
  END SUBROUTINE init_operator_coeffs_cell
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Computes the coefficients that determine the scalar product on the primal grid. This
  !! scalar product depends on the grid geometry only and  is used to formulate the primitive
  !! equations in weak form. The coefficients are applied in module "mo_scalar_product".
  !! The following components of the data type "ocean_patch" are filled:
  !!   edge2cell_coeff  : coefficients for edge to cell mapping
  !!   edge2cell_coeff_t: coefficients for transposed of edge to cell mappings
  !!   edge2vert_coeff  : coefficients for edge to vertex mapping
  !!   edge2vert_coeff_t: coefficients for transposed of edge to vertex mappings
  !!   fixed_vol_norm   : summed volume weight of moved cell
  !!   variable_vol_norm: volume weight at the edges of moved cell
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE init_operator_coeffs_vertex( patch_2D, operators_coefficients)
    TYPE(t_patch), TARGET, INTENT(INOUT) :: patch_2D
    TYPE(t_operator_coeff),INTENT(INOUT) :: operators_coefficients

!Local variables
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc     (1:nproma,1:patch_2D%nblks_v,1:no_dual_edges)
    TYPE(t_cartesian_coordinates) :: edge2vert_coeff_cc_t   (1:nproma,1:patch_2D%nblks_e,1:2)
    TYPE(t_cartesian_coordinates) :: vertex_position, edge_center, vertex_center
    TYPE(t_cartesian_coordinates) :: dist_vector,rot_dist_vector, dist_vector_basic

    REAL(wp)                      :: edge2edge_viavert_coeff(1:nproma,1:patch_2D%nblks_e,1:2*no_dual_edges )
    REAL(wp)                      :: norm, orientation, length, dist_vector_orientedLength
    TYPE(t_cartesian_coordinates) :: z

    INTEGER :: ictr,edge_block_vertex, edge_index_vertex
    INTEGER :: vert_edge
    INTEGER :: edge_block, edge_index
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor
    INTEGER :: level

    TYPE(t_subset_range), POINTER :: owned_edges, all_edges        
    TYPE(t_subset_range), POINTER :: owned_verts         
    !-----------------------------------------------------------------------
    owned_edges => patch_2D%edges%owned
    all_edges   => patch_2D%edges%all
    owned_verts => patch_2D%verts%owned

    edge2vert_coeff_cc(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc(:,:,:)%x(3) = 0.0_wp

    edge2vert_coeff_cc_t(:,:,:)%x(1) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(2) = 0.0_wp
    edge2vert_coeff_cc_t(:,:,:)%x(3) = 0.0_wp
    edge2edge_viavert_coeff(:,:,:)   = 0.0_wp
   !-------------------------------------------
    ! 1) compute:
    !   edge2vert_coeff_cc
    !   variable_dual_vol_norm will be handled in boundary_coeff-sbr
    !variable_dual_vol_norm(:,:,:)  = 0.0_wp

    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index

        vertex_position%x = patch_2D%verts%cartesian(vertex_index, vertex_block)%x

        DO neigbor=1, patch_2D%verts%num_edges(vertex_index,vertex_block)  !no_dual_edges

          edge_index = patch_2D%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch_2D%verts%edge_blk(vertex_index, vertex_block, neigbor)

          IF (edge_block > 0) THEN
            ! we got an adjacent edge
            dist_vector = distance_vector(  &
              & dual_edge_middle(edge_index, edge_block), &
              & vertex_position, &
              & patch_2D%geometry_info)              

            ! the dist_vector has cartesian length
            ! if we use spherical distance we need to recalculate
            ! its length
            IF (LARC_LENGTH) THEN
              ! this will NOT work with the plane geometries !
              length = arc_length(vertex_position, dual_edge_middle(edge_index, edge_block))
              norm = SQRT(SUM( dist_vector%x * dist_vector%x ))
              dist_vector%x = dist_vector%x * length / norm
            ELSE
              length = SQRT(SUM( dist_vector%x * dist_vector%x ))
            ENDIF

            z = get_surface_normal(dual_edge_middle(edge_index, edge_block), patch_2D%geometry_info)
            rot_dist_vector = vector_product(z, dist_vector)
            orientation = DOT_PRODUCT( rot_dist_vector%x,                         &
              & patch_2D%edges%primal_cart_normal(edge_index, edge_block)%x)

            IF ( orientation * patch_2D%verts%edge_orientation(vertex_index, vertex_block, neigbor) <=0) THEN
               write(0,*) "vertex, edgeInVertexList=", vertex_index, vertex_block, neigbor
               write(0,*) "Edge center:", dual_edge_middle(edge_index, edge_block)%x
               write(0,*) "vertex location:", vertex_position%x
               write(0,*) "dist_vector:", dist_vector%x
               write(0,*) "rot_dist_vector:", rot_dist_vector%x
               write(0,*) "orientations", orientation, patch_2D%verts%edge_orientation(vertex_index, vertex_block, neigbor)               
               CALL finish("init_operator_coeffs_vertex", "wrong orientation foredge2vert_coeff_cc" )
            ENDIF
            
            IF (orientation < 0.0_wp) rot_dist_vector%x = - rot_dist_vector%x

            edge2vert_coeff_cc(vertex_index, vertex_block, neigbor)%x = &
              & rot_dist_vector%x *                                     &
              & dual_edge_length(edge_index, edge_block)
              
          ENDIF !(edge_block > 0) THEN
        ENDDO !neigbor=1,6
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch_2D, edge2vert_coeff_cc(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,6
    ! edge2vert_coeff_cc is computed

    !copy 2D to 3D structure
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      DO level = 1, n_zlev
        DO neigbor=1,no_dual_edges

          operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(1)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(1)

          operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(2)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(2)

          operators_coefficients%edge2vert_coeff_cc(:,level,vertex_block,neigbor)%x(3)  &
            &= edge2vert_coeff_cc(:,vertex_block,neigbor)%x(3)

!          operators_coefficients%variable_dual_vol_norm(:,level,vertex_block,neigbor)&
!          &=variable_dual_vol_norm(:,vertex_block,neigbor)
        ENDDO ! neigbor=1,no_dual_edges
      ENDDO  !  level = 1, n_zlev
    ENDDO ! vertex_block
    DO neigbor=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,6
    !----------------------------------------------------

    !----------------------------------------------------
    ! 7) compute:
    !   edge2vert_coeff_cc_t
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = dual_edge_middle(edge_index, edge_block)%x

        DO neigbor=1,2

          edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x = 0.0_wp

          vertex_index = patch_2D%edges%vertex_idx(edge_index, edge_block, neigbor)
          vertex_block = patch_2D%edges%vertex_blk(edge_index, edge_block, neigbor)

          dist_vector = distance_vector(edge_center,  &
            & patch_2D%verts%cartesian(vertex_index, vertex_block), &
            & patch_2D%geometry_info)
            
          edge2vert_coeff_cc_t(edge_index, edge_block, neigbor)%x =           &
            & dist_vector%x *                                                 &
            & patch_2D%edges%tangent_orientation(edge_index, edge_block)  /   &
            & prime_edge_length(edge_index, edge_block)

        ENDDO !neigbor=1,2

      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block
    !-------------------
    ! sync the results
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2
    ! edge2vert_coeff_cc_t is computed

    DO edge_block = owned_edges%start_block, owned_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(owned_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          operators_coefficients%edge2vert_coeff_cc_t(edge_index,level,edge_block,1)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,1)%x

          operators_coefficients%edge2vert_coeff_cc_t(edge_index,level,edge_block,2)%x &
          &=  edge2vert_coeff_cc_t(edge_index,edge_block,2)%x

        ENDDO
      ENDDO
    ENDDO
    DO neigbor=1,2
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(1))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(2))
      CALL sync_patch_array(SYNC_E, patch_2D, edge2vert_coeff_cc_t(:,:,neigbor)%x(3))
    ENDDO ! neigbor=1,2

    !----------------------------------------------------
    ! 8) compute
    !   edge2edge_viavert calculation
    ! For the vorticity advection coeffecients see
    ! B. Perot, Conservation Properties of Unstructured Staggered Mesh Schemes
    !  Eq. (14) and (18)
    DO edge_block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index

        edge_center%x = dual_edge_middle(edge_index, edge_block)%x

        DO neigbor=1,2

          vertex_index   = patch_2D%edges%vertex_idx(edge_index, edge_block, neigbor)
          vertex_block   = patch_2D%edges%vertex_blk(edge_index, edge_block, neigbor)
          vertex_center%x= patch_2D%verts%cartesian(vertex_index, vertex_block)%x

          ! note on orientation: by default the "positive" orientation used in the eqs is:
          !   outwards of the vertex for the tangent, and counterclockwise for the normal
          !   ie. the tangent-normal forms a right-hand system
          !   Note this is ALWAYS in reference to the given vertex
          !       this is the case in our coridinate system when verts%edge_orientation = 1
          !       or, for vertex 1 of the edge, when edges%tangent_orientation = 1  
          !dist_vector_basic%x = (edge_center%x - vertex_center%x) &

          dist_vector_basic = distance_vector(vertex_center, edge_center, patch_2D%geometry_info)
          dist_vector_basic%x = dist_vector_basic%x & ! opposite orientation since we use left-hand system
            & * (3 - 2 * neigbor) * patch_2D%edges%tangent_orientation(edge_index, edge_block) &
!            & * dual_edge_length(edge_index, edge_block)          & not used in a not-flux form
            &  / prime_edge_length(edge_index, edge_block)
          ! this is multiplied by the dual edge length on unit sphere

          !IF(neigbor==1)ictr = 0
          !IF(neigbor==2)ictr = no_dual_edges

          ictr = (neigbor - 1)*no_dual_edges

!           DO vert_edge=1,patch_2D%verts%num_edges(vertex_index,vertex_block)!no_dual_edges
!             ictr=ictr+1
!             !actual edge
!             edge_index_cell = patch_2D%verts%edge_idx(vertex_index, vertex_block, vert_edge)
!             edge_block_cell = patch_2D%verts%edge_blk(vertex_index, vertex_block, vert_edge)
!             dist_vector%x  =  dual_edge_middle(edge_index_cell, edge_block_cell)%x - vertex_center%x
! 
!             dist_vector = vector_product(dist_vector, dual_edge_middle(edge_index_cell, edge_block_cell))
!             orientation = DOT_PRODUCT( dist_vector%x,                         &
!                & patch_2D%edges%primal_cart_normal(edge_index_cell, edge_block_cell)%x)
!             ! orientation should not be 0, since this would mean that the prime and dual are parallel
!             ! overall this calculation should be derived from the verts%edge_orientation
!             ! orientation will receive a value -1, or 1 based on the previous,
!             ! then multuplied by -1 if neigbor=2, otherwise unchanged
!             orientation = SIGN(1.0_wp,orientation) * (3.0_wp - 2.0_wp * REAL(neigbor,wp))
! 
!             !The dot product is the cosine of the angle between vectors from dual cell centers
!             !to dual cell edges 
!             edge2edge_viavert_coeff(edge_index,edge_block,ictr)         &
!               & = orientation                                           &
!               & * DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)        &
!               & * patch_2D%edges%tangent_orientation(edge_index, edge_block)&
!               & * (dual_edge_length(edge_index_cell, edge_block_cell)   &
!               &    / prime_edge_length(edge_index, edge_block))
! 
!           END DO
          
          DO vert_edge=1,patch_2D%verts%num_edges(vertex_index,vertex_block)!no_dual_edges
            ictr=ictr+1
            !actual edge
            edge_index_vertex = patch_2D%verts%edge_idx(vertex_index, vertex_block, vert_edge)
            edge_block_vertex = patch_2D%verts%edge_blk(vertex_index, vertex_block, vert_edge)
            
            IF (edge_index == edge_index_vertex .and. edge_block == edge_block_vertex) THEN
              ! the result is 0, since the external product (see below) of this edge is
              ! perpedicular to itself, and the dot product is 0
              edge2edge_viavert_coeff(edge_index,edge_block,ictr) = 0.0_wp
            ELSE

              dist_vector = distance_vector(  &
                & dual_edge_middle(edge_index_vertex, edge_block_vertex), &
                & vertex_center, patch_2D%geometry_info)

              dist_vector%x  = dist_vector%x  &
                & * patch_2D%verts%edge_orientation(vertex_index, vertex_block, vert_edge)

              ! we need the normalized half of the dual_edge_middle + vertex at z
              z = planar_middle( &
               & dual_edge_middle(edge_index_vertex, edge_block_vertex), &
               & vertex_center, patch_2D%geometry_info)
              z = get_surface_normal(z, patch_2D%geometry_info) ! get_surface_normal has to do z%x=z%x/SQRT(SUM(z%x * z%x))
              z%x=z%x/SQRT(SUM(z%x * z%x))
              dist_vector = vector_product(dist_vector, z)
              ! the dist_vector has still dual_edge_middle-vertex length

!               dist_vector%x = &
!                 & - patch_2D%edges%primal_cart_normal(edge_index_vertex, edge_block_vertex)%x &
!                 & * distance(  &
!                 &     dual_edge_middle(edge_index_vertex, edge_block_vertex), &
!                 &     vertex_center, patch_2D%geometry_info) &
!                 & * patch_2D%verts%edge_orientation(vertex_index, vertex_block, vert_edge)

              ! adjust the orientation along the vn of the edge
              edge2edge_viavert_coeff(edge_index,edge_block,ictr)             &
                & =  DOT_PRODUCT(dist_vector_basic%x,dist_vector%x)           &
                &  * dual_edge_length(edge_index_vertex, edge_block_vertex)

! 
!               write(0,*) dist_vector%x
!               write(0,*) patch_2D%edges%primal_cart_normal(edge_index_vertex, edge_block_vertex)%x * &
!                   dist_vector_orientedLength
!               write(0,*) edge2edge_viavert_coeff(edge_index,edge_block,ictr),  &
!                 DOT_PRODUCT(dist_vector_basic%x, patch_2D%edges%primal_cart_normal(edge_index_vertex, edge_block_vertex)%x) * &
!                 SQRT(SUM(dist_vector%x * dist_vector%x)) *  patch_2D%verts%edge_orientation(vertex_index, vertex_block, vert_edge) &
!                   &  * dual_edge_length(edge_index_vertex, edge_block_vertex)
            ENDIF

          END DO
        ENDDO !neigbor=1,2
      ENDDO ! edge_index = start_index, end_index
    ENDDO ! edge_block = owned_edges%start_block, owned_edges%end_block

    !-------------------
    DO ictr=1, 2*no_dual_edges
!      write(0,*)'ictr:',ictr
      CALL sync_patch_array(SYNC_E, patch_2D, edge2edge_viavert_coeff(:,:,ictr))
    ENDDO

    DO edge_block = all_edges%start_block, all_edges%end_block
      DO level = 1, n_zlev
        CALL get_index_range(all_edges, edge_block, start_index, end_index)
        DO edge_index =  start_index, end_index

          operators_coefficients%edge2edge_viavert_coeff(edge_index,level,edge_block,1:2*no_dual_edges) = &
            & edge2edge_viavert_coeff(edge_index,edge_block,1:2*no_dual_edges)

        ENDDO
      ENDDO
    ENDDO
!    DO neigbor=1,2*no_dual_edges
!      CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2edge_viavert_coeff(:,:,:,neigbor))
!    END DO
   !-------------------------------------------
  END SUBROUTINE init_operator_coeffs_vertex
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !> Modify coefficients at the boundary
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE apply_boundary2coeffs( patch_3D, operators_coefficients)
    ! !
    TYPE(t_patch_3D ),TARGET, INTENT(INOUT) :: patch_3D
    TYPE(t_operator_coeff), INTENT(INOUT)   :: operators_coefficients

    !Local variables
    INTEGER :: jk, jc, block, je, ibe, ile, jev, jv,ie
    INTEGER :: il_e, ib_e, il_c, ib_c, ictr
    INTEGER :: edges_startidx, edges_endidx
    INTEGER :: cells_startidx, cells_endidx
    INTEGER :: i_startidx_v, i_endidx_v

    INTEGER :: sea_edges_per_vertex(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: ivertex_bnd_edge_idx(4), ivertex_bnd_edge_blk(4)  !maximal 4 boundary edges in a dual loop.
    INTEGER :: i_edge_idx(4)
    REAL(wp) :: zarea_fraction(nproma,n_zlev,patch_3D%p_patch_2D(1)%nblks_v)
    INTEGER :: cell_idx_1, cell_blk_1
    INTEGER :: cell_idx_2, cell_blk_2
    INTEGER :: boundary_counter
    INTEGER :: neigbor, k, boundary_level, i, max_level
    INTEGER :: cell_index, cell_block, edge_index_of_cell, edge_block_of_cell, k_coeff

    !INTEGER :: vertex_edge
    TYPE(t_cartesian_coordinates) :: cell1_cc, cell2_cc, vertex_cc
    REAL(wp) :: grid_radius_squared
    INTEGER, PARAMETER :: MIN_SEA_BOUNDARYLEVEL = -99999
    TYPE(t_subset_range), POINTER :: all_edges, owned_edges
    TYPE(t_subset_range), POINTER :: all_cells, owned_cells
    TYPE(t_subset_range), POINTER :: owned_verts!, in_domain_verts
    TYPE(t_patch), POINTER        :: patch_2D
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:apply_boundary2coeffs')
    !-----------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    patch_2D    => patch_3D%p_patch_2D(1)
    all_cells   => patch_2D%cells%all
    all_edges   => patch_2D%edges%all
    owned_cells => patch_2D%cells%owned
    owned_edges => patch_2D%edges%owned
    owned_verts => patch_2D%verts%owned
    !in_domain_verts  => patch_3D%p_patch_2D(1)%verts%in_domain

    sea_edges_per_vertex(:,:,:)                       = 0
    grid_radius_squared = grid_sphere_radius * grid_sphere_radius

    ! calculate cells SeaBoundaryLevel
    ! first initialize levels 1,-1
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, cells_startidx, cells_endidx)
      DO jc = cells_startidx, cells_endidx
        DO jk = 1, n_zlev
          IF (patch_3D%lsm_c(jc,jk,block) > boundary) THEN
            operators_coefficients%cells_SeaBoundaryLevel(jc,jk,block) = 1 ! land
          ELSEIF ( patch_3D%lsm_c(jc,jk,block) == sea_boundary ) THEN
            operators_coefficients%cells_SeaBoundaryLevel(jc,jk,block) = -1 ! sea boundary
          ELSE
            operators_coefficients%cells_SeaBoundaryLevel(jc,jk,block) = MIN_SEA_BOUNDARYLEVEL ! inner sea
          ENDIF
        ENDDO
      END DO
    END DO

    ! calclulate the next levels
    DO boundary_level=1,2 ! max calculated level is -3

      DO block = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, block, cells_startidx, cells_endidx)
        DO jc = cells_startidx, cells_endidx
          DO jk = 1, n_zlev
            IF (operators_coefficients%cells_SeaBoundaryLevel(jc,jk,block) == MIN_SEA_BOUNDARYLEVEL) THEN
              ! find the max level of the neigbors
              max_level = MIN_SEA_BOUNDARYLEVEL
              DO i=1,patch_2D%cells%num_edges(jc,block)
                cell_idx_1 = patch_2D%cells%neighbor_idx(jc,block,i)
                cell_blk_1 = patch_2D%cells%neighbor_blk(jc,block,i)
                max_level = MAX(max_level,operators_coefficients%cells_SeaBoundaryLevel(cell_idx_1,jk,cell_blk_1))
              ENDDO
             
              IF (max_level /= MIN_SEA_BOUNDARYLEVEL) THEN
                ! if we are in deep sea do nothing, else compute the level
                IF (max_level > -boundary_level) THEN ! we have a problem
                  write(0,*) -boundary_level,  max_level
                  CALL finish("calculate cells SeaBoundaryLevel", "max_level > -boundary_level")
                ENDIF
                operators_coefficients%cells_SeaBoundaryLevel(jc,jk,block) = max_level - 1
              ENDIF

            ENDIF

          ENDDO
        END DO
      END DO

      CALL sync_patch_array(SYNC_C, patch_2D, operators_coefficients%cells_SeaBoundaryLevel)

    ENDDO

    ! calculate edges_SeaBoundaryLevel based on cells_SeaBoundaryLevel
    ! first calculate levels 1,0,-1
    DO block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, block, edges_startidx, edges_endidx)
      DO je = edges_startidx, edges_endidx
        cell_idx_1 = patch_2D%edges%cell_idx(je,block,1)
        cell_blk_1 = patch_2D%edges%cell_blk(je,block,1)
        cell_idx_2 = patch_2D%edges%cell_idx(je,block,2)
        cell_blk_2 = patch_2D%edges%cell_blk(je,block,2)

        DO jk = 1, n_zlev
          IF (patch_3D%lsm_e(je,jk,block) > boundary) THEN
            operators_coefficients%edges_SeaBoundaryLevel(je,jk,block) = 1 ! land
          ELSEIF ( patch_3D%lsm_e(je,jk,block) == boundary ) THEN
            operators_coefficients%edges_SeaBoundaryLevel(je,jk,block) = 0 ! boundary
          ELSE
            operators_coefficients%edges_SeaBoundaryLevel(je,jk,block) = &
              & MAX(operators_coefficients%cells_SeaBoundaryLevel(cell_idx_1, jk, cell_blk_1), &
              &     operators_coefficients%cells_SeaBoundaryLevel(cell_idx_2, jk, cell_blk_2))
            IF (operators_coefficients%edges_SeaBoundaryLevel(je,jk,block) > -1) &
              & CALL finish("edges_SeaBoundaryLevel","edges_SeaBoundaryLevel > -1")
          ENDIF
        ENDDO
      END DO
    END DO
    CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edges_SeaBoundaryLevel)

    !-------------------------------------------------------------
    !0. check the coefficients for edges, these are:
    !     grad_coeff
    !     edge2cell_coeff_cc_t
    !     edge2vert_coeff_cc_t
    !
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, edges_startidx, edges_endidx)
      DO je = edges_startidx, edges_endidx
        DO jk = 1, n_zlev

          IF ( patch_3D%lsm_e(je,jk,block) /= sea ) THEN
            operators_coefficients%grad_coeff          (je,jk,block) = 0.0_wp
            operators_coefficients%edge2cell_coeff_cc_t(je,jk,block,1)%x(1:3) = 0.0_wp
            operators_coefficients%edge2cell_coeff_cc_t(je,jk,block,2)%x(1:3) = 0.0_wp
            operators_coefficients%edge2vert_coeff_cc_t(je,jk,block,1)%x(1:3) = 0.0_wp
            operators_coefficients%edge2vert_coeff_cc_t(je,jk,block,2)%x(1:3) = 0.0_wp
          ENDIF
        ENDDO
      END DO
    END DO
    !-------------------------------------------------------------
    !All the coefficients "edge2edge_viacell_coeff" are set to zero
    !if the actual edge under consideration is an land edge.
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, edges_startidx, edges_endidx)
      DO jk = 1, n_zlev
        DO je = edges_startidx, edges_endidx
          IF ( patch_3D%lsm_e(je,jk,block) /= sea ) THEN
            operators_coefficients%edge2edge_viacell_coeff(je,jk,block,:) = 0.0_wp
            operators_coefficients%edge2edge_viavert_coeff(je,jk,block,:) = 0.0_wp
          ENDIF
        END DO
      END DO
    END DO
    !-------------------------------------------------------------
    !The coefficients "edge2edge_viacell_coeff" are set to zero
    !for which the stencil contains is an land edge.
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, edges_startidx, edges_endidx)
      DO jk = 1, n_zlev
        DO je = edges_startidx, edges_endidx

          !Handle neighbour cell 1
          ictr  = 0
          il_c  = patch_2D%edges%cell_idx(je,block,1)
          ib_c  = patch_2D%edges%cell_blk(je,block,1)
          IF (il_c > 0) THEN
            DO ie=1, no_primal_edges
              ictr =ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)
              IF (il_e > 0) THEN
                IF ( patch_3D%lsm_e(il_e,jk,ib_e) /= sea ) THEN
                  operators_coefficients%edge2edge_viacell_coeff(je,jk,block,ictr)=0.0_wp
                ENDIF
              ENDIF
            END DO
          ENDIF

          !Handle neighbour cell 2
          ictr  = no_primal_edges
          il_c  = patch_2D%edges%cell_idx(je,block,2)
          ib_c  = patch_2D%edges%cell_blk(je,block,2)
          IF (il_c > 0) THEN
            DO ie=1, no_primal_edges
              ictr =ictr+1
              il_e = patch_2D%cells%edge_idx(il_c,ib_c,ie)
              ib_e = patch_2D%cells%edge_blk(il_c,ib_c,ie)

              IF (il_e > 0) THEN
                IF ( patch_3D%lsm_e(il_e,jk,ib_e) /= sea ) THEN
                  operators_coefficients%edge2edge_viacell_coeff(je,jk,block,ictr)=0.0_wp
                ENDIF
              ENDIF
            END DO
          ENDIF
        END DO
      END DO
    END DO 

    !-------------------------------------------------------------
    ! Normalize "edge2edge_viacell_coeff" 
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, edges_startidx, edges_endidx)
      DO jk = 1, n_zlev
        DO je = edges_startidx, edges_endidx
          IF ( patch_3D%lsm_e(je,jk,block) == sea ) THEN
            cell_idx_1 = patch_2D%edges%cell_idx(je,block,1)
            cell_blk_1 = patch_2D%edges%cell_blk(je,block,1)

            cell_idx_2 = patch_2D%edges%cell_idx(je,block,2)
            cell_blk_2 = patch_2D%edges%cell_blk(je,block,2)

            operators_coefficients%edge2edge_viacell_coeff(je,jk,block,1:no_primal_edges) &
            &= operators_coefficients%edge2edge_viacell_coeff(je,jk,block,1:no_primal_edges)&
            &/operators_coefficients%fixed_vol_norm(cell_idx_1,jk,cell_blk_1)

            operators_coefficients%edge2edge_viacell_coeff(je,jk,block,no_primal_edges+1:2*no_primal_edges) &
            &= operators_coefficients%edge2edge_viacell_coeff(je,jk,block,no_primal_edges+1:2*no_primal_edges)&
            &/operators_coefficients%fixed_vol_norm(cell_idx_2,jk,cell_blk_2)
          ENDIF

        END DO
      END DO
    END DO 

    !-------------------------------------------------------------
    !Fill edge2edge_viacell_coeff_top and edge2edge_viacell_coeff_integrated from edge2edge_viacell_coeff
    !
    ! the top is just the first level rearranged:
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, edges_startidx, edges_endidx)
      DO je = edges_startidx, edges_endidx
        DO k=1, 2 * no_primal_edges
          operators_coefficients%edge2edge_viacell_coeff_top(k, je, block) = &
             operators_coefficients%edge2edge_viacell_coeff(je, 1, block, k)
        ENDDO
!         write(0,*) patch_3D%lsm_e(je,1,block), ", edge2edge_viacell:", operators_coefficients%edge2edge_viacell_coeff_top(:, je, block)
      ENDDO
    ENDDO

    ! the integrated are the rest of levels > 1, weighted by prism_thick_e and integrated:
    DO block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, block, edges_startidx, edges_endidx)
      DO je = edges_startidx, edges_endidx

        ! zero for two cells/ six edges
        DO k_coeff = 1, 2 * no_primal_edges
           operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, block) = 0.0_wp
        ENDDO

        ! the first cell
        cell_index   = patch_2D%edges%cell_idx(je,block,1)
        cell_block = patch_2D%edges%cell_blk(je,block,1)
        IF (cell_index > 0) THEN ! in case it's a lateral boundary edge
          DO k=1, no_primal_edges
            edge_index_of_cell = patch_2D%cells%edge_idx(cell_index, cell_block, k)
            edge_block_of_cell = patch_2D%cells%edge_blk(cell_index, cell_block, k)
            k_coeff = k
            DO jk=2, patch_3d%p_patch_1d(1)%dolic_e(je,block)
              operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, block) = &
                 operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, block) + &
                 operators_coefficients%edge2edge_viacell_coeff(je, jk, block, k_coeff) * &
                 patch_3D%p_patch_1D(1)%prism_thick_e(edge_index_of_cell, jk, edge_block_of_cell)

           !   write(0,*) block, je, jk, k_coeff, operators_coefficients%edge2edge_viacell_coeff(je, jk, block, k_coeff), &
           !     & patch_3D%p_patch_1D(1)%prism_thick_e(edge_index_of_cell, jk, edge_block_of_cell), &
           !    &  operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, block)
            ENDDO
          ENDDO
        ENDIF

        ! the second cell
        cell_index   = patch_2D%edges%cell_idx(je,block,2)
        cell_block = patch_2D%edges%cell_blk(je,block,2)
        IF (cell_index > 0) THEN ! in case it's a lateral boundary edge
          DO k=1, no_primal_edges
            edge_index_of_cell = patch_2D%cells%edge_idx(cell_index, cell_block, k)
            edge_block_of_cell = patch_2D%cells%edge_blk(cell_index, cell_block, k)
            k_coeff = no_primal_edges + k
            DO jk=2, patch_3d%p_patch_1d(1)%dolic_e(je,block)
              operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, block) = &
                 operators_coefficients%edge2edge_viacell_coeff_integrated(k_coeff, je, block) + &
                 operators_coefficients%edge2edge_viacell_coeff(je, jk, block, k_coeff) * &
                 patch_3D%p_patch_1D(1)%prism_thick_e(edge_index_of_cell, jk, edge_block_of_cell)

            ENDDO
          ENDDO
        ENDIF

      ENDDO
    ENDDO

    ! zeros edge2edge_viacell_coeff_all
!    DO block = all_edges%start_block, all_edges%end_block
!      CALL get_index_range(all_edges, block, edges_startidx, edges_endidx)
!      DO je = edges_startidx, edges_endidx
!        DO k=1, 2 * no_primal_edges
!          operators_coefficients%edge2edge_viacell_coeff_all(k, je, block) = 0.0_wp
!        ENDDO
!      ENDDO
!    ENDDO
    operators_coefficients%edge2edge_viacell_coeff_all(:,:,:) = 0.0_wp
    !-------------------------------------------------------------
    !1) Set coefficients for div and grad to zero at boundary edges
    ! Also for operators_coefficients%edge2cell_coeff_cc
    DO block = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, block, cells_startidx, cells_endidx)
      DO jk=1,n_zlev
        DO jc = cells_startidx, cells_endidx
          DO je = 1, patch_2D%cells%num_edges(jc,block)!no_primal_edges!

            ile = patch_2D%cells%edge_idx(jc,block,je)
            ibe = patch_2D%cells%edge_blk(jc,block,je)

            IF ( patch_3D%lsm_e(ile,jk,ibe) /= sea ) THEN
              operators_coefficients%div_coeff(jc,jk,block,je) = 0.0_wp
              operators_coefficients%edge2cell_coeff_cc(jc,jk,block,je)%x(1:3) = 0.0_wp
            ENDIF
          ENDDO ! je = 1, patch_2D%cells%num_edges(jc,block)
        ENDDO ! jc = cells_startidx, cells_endidx
      END DO ! jk=1,n_zlev
    END DO ! block = all_cells%start_block, all_cells%end_block
    !-------------------------------------------------------------
    !2) prepare coefficients for rot at boundary edges
    ! this is done on owned vertices, so sync is required
    DO block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, block, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v

          ivertex_bnd_edge_idx(1:4)      = 0
          ivertex_bnd_edge_blk(1:4)      = 0
          !i_edge_idx(1:4)         = 0
          !z_orientation(1:4)      = 0.0_wp
          boundary_counter        = 0

          DO jev = 1, patch_2D%verts%num_edges(jv,block)

            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,block,jev)
            ibe = patch_2D%verts%edge_blk(jv,block,jev)
            !Check, if edge is sea or boundary edge and take care of dummy edge
            ! edge with indices ile, ibe is sea edge
            ! edge with indices ile, ibe is boundary edge

            IF ( patch_3D%lsm_e(ile,jk,ibe) == SEA ) THEN
              sea_edges_per_vertex(jv,jk,block) = sea_edges_per_vertex(jv,jk,block) + 1
            ELSEIF ( patch_3D%lsm_e(ile,jk,ibe) == BOUNDARY ) THEN

              !increase boundary edge counter
              boundary_counter = boundary_counter + 1

              operators_coefficients%bnd_edges_per_vertex(jv,jk,block) &
                & = operators_coefficients%bnd_edges_per_vertex(jv,jk,block) +1

              IF (boundary_counter > 4) THEN
                !maximal 4 boundary edges per dual loop are allowed: somethings wrong with the grid
                CALL message (TRIM('sbr nonlinear Coriolis'), &
                  & 'more than 4 boundary edges per dual loop: something is wrong with the grid')
                CALL finish (routine,'Grid-boundary error !!')
              ENDIF
              ivertex_bnd_edge_idx(boundary_counter) = ile
              ivertex_bnd_edge_blk(boundary_counter) = ibe
              !z_orientation(boundary_counter) = patch_2D%verts%edge_orientation(jv,block,jev)
              i_edge_idx(boundary_counter)    = jev

              operators_coefficients%vertex_bnd_edge_idx(jv,jk,block,boundary_counter)= ile
              operators_coefficients%vertex_bnd_edge_blk(jv,jk,block,boundary_counter)= ibe
!               operators_coefficients%orientation(jv,jk,block,boundary_counter) = &
!                 & patch_2D%verts%edge_orientation(jv,block,jev)
              operators_coefficients%boundaryEdge_Coefficient_Index(jv,jk,block,boundary_counter)    = jev

            END IF
          END DO ! jev = 1, patch_2D%verts%num_edges(jv,block)

          IF( MOD(boundary_counter,2) /= 0 ) THEN
            CALL finish (routine,'MOD(boundary_counter,2) /= 0 !!')
          ENDIF

          !---------------------------------------------------------------------------------
          !Modified area calculation
          vertex_cc = patch_2D%verts%cartesian(jv,block)
          DO jev = 1, patch_2D%verts%num_edges(jv,block)
            ! get line and block indices of edge jev around vertex jv
            ile = patch_2D%verts%edge_idx(jv,block,jev)
            ibe = patch_2D%verts%edge_blk(jv,block,jev)
            !get neighbor cells
            cell_idx_1 = patch_2D%edges%cell_idx(ile,ibe,1)
            cell_idx_2 = patch_2D%edges%cell_idx(ile,ibe,2)
            cell_blk_1 = patch_2D%edges%cell_blk(ile,ibe,1)
            cell_blk_2 = patch_2D%edges%cell_blk(ile,ibe,2)

!             IF ( patch_3D%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
!               cell1_cc%x  = patch_2D%cells%cartesian_center(cell_idx_1,cell_blk_1)%x
!               cell2_cc%x  = patch_2D%cells%cartesian_center(cell_idx_2,cell_blk_2)%x
! 
!               !Check, if edge is sea or boundary edge and take care of dummy edge
!               !edge with indices ile, ibe is sea edge
!               !Add up for wet dual area.
!               !IF ( v_base%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
!               operators_coefficients%variable_dual_vol_norm(jv,jk,block,jev)= &
!                 & planar_triangle_area(cell1_cc, vertex_cc, cell2_cc, patch_2D%geometry_info) ! this need to be on the plane area for the sphere
!               ! edge with indices ile, ibe is boundary edge
!             ELSE IF ( patch_3D%lsm_e(ile,jk,ibe) == boundary ) THEN
!               operators_coefficients%variable_dual_vol_norm(jv,jk,block,jev)=0.0_wp
!               !0.5_wp*planar_triangle_area(cell1_cc, vertex_cc, cell2_cc)
!             END IF
          END DO

          !---------------------------------------------------------------------------------------------
!           DO je = 1, boundary_counter
!             ivertex_bnd_edge_idx(je) = operators_coefficients%vertex_bnd_edge_idx(jv,jk,block,je)
!             ivertex_bnd_edge_blk(je) = operators_coefficients%vertex_bnd_edge_blk(jv,jk,block,je)
! 
!             ! needs to be re-examined !
!             operators_coefficients%rot_coeff(jv,jk,block,i_edge_idx(je) )=&
!               & 0.5_wp*patch_2D%edges%tangent_orientation(ivertex_bnd_edge_idx(je),ivertex_bnd_edge_blk(je)) * &
!               & prime_edge_length(ivertex_bnd_edge_idx(je),ivertex_bnd_edge_blk(je)) * grid_sphere_radius
!               ! this is the real distance on the Earth
! 
!           ENDDO
        END DO ! jv = i_startidx_v, i_endidx_v

      END DO ! jk = 1, n_zlev
    END DO ! block = owned_verts%start_block, owned_verts%end_block
    ! sync rot_coeff is done after area normalization at the end
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !The dynamical changing coefficient for the surface layer
!     DO block = all_cells%start_block, all_cells%end_block
!       CALL get_index_range(all_cells, block, cells_startidx, cells_endidx)
!       DO jc = cells_startidx, cells_endidx
! 
!         DO je = 1, patch_2D%cells%num_edges(jc,block)
!           operators_coefficients%edge2cell_coeff_cc_dyn(jc,1,block,je)%x = &
!             operators_coefficients%edge2cell_coeff_cc(jc,1,block,je)%x
!         ENDDO 
!       END DO ! jc = cells_startidx, cells_endidx
!     END DO ! block = all_cells%start_block, all_cells%end_block
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !3.3) Edge to vert coefficient
    DO block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, block, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v

          DO je = 1, patch_2D%verts%num_edges(jv,block)
            ile = patch_2D%verts%edge_idx(jv,block,je)
            ibe = patch_2D%verts%edge_blk(jv,block,je)

             IF ( patch_3D%lsm_e(ile,jk,ibe) /= sea) THEN
              operators_coefficients%edge2vert_coeff_cc(jv,jk,block,je)%x(1:3) = 0.0_wp
!               operators_coefficients%variable_dual_vol_norm(jv,jk,block,je)    = 0.0_wp
            ENDIF
          ENDDO ! je = 1, patch_2D%verts%num_edges(jv,block)
        ENDDO ! jv = i_startidx_v, i_endidx_v
      END DO ! jk = 1, n_zlev
    END DO ! block = owned_verts%start_block, owned_verts%end_block
    ! sync the result
!     DO je=1,no_dual_edges
!       CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%variable_dual_vol_norm(:,:,:, je))
!     ENDDO
    !-------------------------------------------------------------

    !-------------------------------------------------------------
    !Merge dual area calculation with coefficients
    ! note: sea_edges_per_vertex has been calculated on the owned_verts
    !       it does not need to be synced
    zarea_fraction(:,:,:) = 0.0_wp
    DO block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, block, i_startidx_v, i_endidx_v)
      DO jk = 1, n_zlev
        DO jv = i_startidx_v, i_endidx_v
          zarea_fraction(jv,jk,block) = 0.0_wp

!           IF ( sea_edges_per_vertex(jv,jk,block) == no_dual_edges ) THEN ! we have to count for lateral boundaries at the top
!             zarea_fraction(jv,jk,block)= patch_2D%verts%dual_area(jv,block) / grid_radius_squared
!             !zarea_fraction(jv,jk,block)=SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,block,:))
! 
!             !ELSEIF(operators_coefficients%bnd_edges_per_vertex(jv,jk,block)/=0)THEN!boundary edges are involved
!           ELSEIF ( sea_edges_per_vertex(jv,jk,block) /= 0 ) THEN
          IF ( sea_edges_per_vertex(jv,jk,block) /= 0 ) THEN

            !Modified area calculation
            vertex_cc = patch_2D%verts%cartesian(jv,block)
            DO jev = 1, patch_2D%verts%num_edges(jv,block)
              ! get line and block indices of edge jev around vertex jv
              ile = patch_2D%verts%edge_idx(jv,block,jev)
              ibe = patch_2D%verts%edge_blk(jv,block,jev)
              !get neighbor cells
              cell_idx_1 = patch_2D%edges%cell_idx(ile,ibe,1)
              cell_idx_2 = patch_2D%edges%cell_idx(ile,ibe,2)
              cell_blk_1 = patch_2D%edges%cell_blk(ile,ibe,1)
              cell_blk_2 = patch_2D%edges%cell_blk(ile,ibe,2)

              !Check, if edge is sea or boundary edge and take care of dummy edge
              !edge with indices ile, ibe is sea edge
              !Add up for wet dual area.
              ! Note that this should be modified.
              !   sea_boundary means an open boundary
              !   boundary means that only the sea cell are should be added
              IF ( patch_3D%lsm_e(ile,jk,ibe) <= sea_boundary ) THEN
                cell1_cc%x  = patch_2D%cells%cartesian_center(cell_idx_1,cell_blk_1)%x
                cell2_cc%x  = patch_2D%cells%cartesian_center(cell_idx_2,cell_blk_2)%x
                zarea_fraction(jv,jk,block) = zarea_fraction(jv,jk,block)  &
                  & + planar_triangle_area(cell1_cc, vertex_cc, cell2_cc, patch_2D%geometry_info)
                ! edge with indices ile, ibe is boundary edge                
              ELSE IF ( patch_3D%lsm_e(ile,jk,ibe) == boundary ) THEN
                ! at least one of the two cells exists and is sea cell
                IF (cell_idx_2 <= 0) THEN
                  cell1_cc%x  = patch_2D%cells%cartesian_center(cell_idx_1,cell_blk_1)%x
                ELSE IF (cell_idx_1 <= 0) THEN
                  cell1_cc%x  = patch_2D%cells%cartesian_center(cell_idx_2,cell_blk_2)%x
                ELSE IF (patch_3D%lsm_c(cell_idx_1,jk,cell_blk_1) <= sea_boundary) THEN
                  cell1_cc%x  = patch_2D%cells%cartesian_center(cell_idx_1,cell_blk_1)%x
                ELSE
                  cell1_cc%x  = patch_2D%cells%cartesian_center(cell_idx_2,cell_blk_2)%x
                ENDIF
!                 zarea_fraction(jv,jk,block) = zarea_fraction(jv,jk,block)  &
!                   & + 0.5_wp*planar_triangle_area(cell1_cc, vertex_cc, cell2_cc)
                ! add only the sea dual area, ie the triagle area between
                ! the vertex, edge centre, and sea cell center
                zarea_fraction(jv,jk,block) = zarea_fraction(jv,jk,block)  &
                  & + planar_triangle_area(cell1_cc, vertex_cc, patch_2D%edges%cartesian_center(ile,ibe), patch_2D%geometry_info)
              ENDIF
              
            END DO ! jev = 1, patch_2D%verts%num_edges(jv,block)
            
            ! this is calculated already on the unit sphere
            ! zarea_fraction(jv,jk,block) = zarea_fraction(jv,jk,block) / grid_radius_squared
            
          ENDIF !( sea_edges_per_vertex(jv,jk,block) == patch_2D%verts%num_edges(jv,block) )
          !The two quantities: 
          !zarea_fraction(jv,jk,block)*(earth_radius*earth_radius)
          !and 
          !patch_2D%verts%dual_area(jv,block)
          !are identical

          !Final coefficient calculation
!           IF (operators_coefficients%bnd_edges_per_vertex(jv,jk,block) > 0) THEN
!             DO jev = 1, no_dual_edges
!               operators_coefficients%rot_coeff(jv,jk,block,jev)=0.0_wp
!             END DO

          IF(zarea_fraction(jv,jk,block) > 0.0_wp)THEN
            
            DO jev = 1, no_dual_edges
              operators_coefficients%edge2vert_coeff_cc(jv,jk,block,jev)%x(1:3)&
                & =operators_coefficients%edge2vert_coeff_cc(jv,jk,block,jev)%x(1:3)/zarea_fraction(jv,jk,block)
                !SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,block,:))!
              operators_coefficients%rot_coeff(jv,jk,block,jev)&
                &=operators_coefficients%rot_coeff(jv,jk,block,jev)/(zarea_fraction(jv,jk,block)*grid_radius_squared)

!               IF (ABS(operators_coefficients%rot_coeff(jv,jk,block,jev)) > 1.0E-2_wp) THEN
!                 write(0,*) "rot_coeff > 1.0E-2_wp: ", operators_coefficients%rot_coeff(jv,jk,block,jev), &
!                   zarea_fraction(jv,jk,block), grid_radius_squared
!               ENDIF

            END DO

          ELSE

            DO jev = 1, no_dual_edges
              operators_coefficients%edge2vert_coeff_cc(jv,jk,block,jev)%x(1:3)=0.0_wp
              operators_coefficients%rot_coeff(jv,jk,block,jev)=0.0_wp
            END DO

          ENDIF
         !!ENDIF !( sea_edges_per_vertex(jv,jk,block) == patch_2D%verts%num_edges(jv,block) )

        ENDDO!jv = i_startidx_v, i_endidx_v
      END DO!jk = 1, n_zlev
    END DO!block = owned_verts%start_block, owned_verts%end_block
    ! sync the result
    DO jev=1,no_dual_edges
      DO jk = 1, n_zlev
        CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,jk,:, jev)%x(1))
        CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,jk,:, jev)%x(2))
        CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%edge2vert_coeff_cc(:,jk,:, jev)%x(3))
      ENDDO
    ENDDO
    DO je=1,no_dual_edges
      CALL sync_patch_array(SYNC_V, patch_2D, operators_coefficients%rot_coeff(:,:,:, je))
    ENDDO
!     DO je=1,no_dual_edges
!       CALL dbg_print('rot_coeff'    ,operators_coefficients%rot_coeff(:,:,:, je), this_mod_name,1, &
!             patch_2D%verts%owned )
!     ENDDO

!     DO jev=1,2*no_dual_edges
!       DO jk = 1, n_zlev
!         CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2edge_viavert_coeff(:,jk,:, jev))
!       ENDDO
!     ENDDO
    CALL sync_patch_array(SYNC_V, patch_2D, zarea_fraction(:,:,:))
    
    DO block = owned_edges%start_block, owned_edges%end_block
      CALL get_index_range(owned_edges, block, edges_startidx, edges_endidx)
      DO je = edges_startidx, edges_endidx
        DO jk = 1, n_zlev
          
          IF ( patch_3D%lsm_e(je,jk,block) <= sea_boundary ) THEN

            DO neigbor=1,2

              jv  = patch_2D%edges%vertex_idx(je, block, neigbor)
              jev = patch_2D%edges%vertex_blk(je, block, neigbor)

              IF(neigbor==1)THEN
                operators_coefficients%edge2edge_viavert_coeff(je,jk,block,1:no_dual_edges)&
                &=operators_coefficients%edge2edge_viavert_coeff(je,jk,block,1:no_dual_edges)&
                &/zarea_fraction(jv,jk,jev)!SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,jev,:))
              ELSEIF(neigbor==2)THEN
                operators_coefficients%edge2edge_viavert_coeff(je,jk,block,no_dual_edges+1:2*no_dual_edges)&
                &=operators_coefficients%edge2edge_viavert_coeff(je,jk,block,no_dual_edges+1:2*no_dual_edges)&
                &/zarea_fraction(jv,jk,jev)!SUM(operators_coefficients%variable_dual_vol_norm(jv,jk,jev,:))
              ENDIF
            
            END DO !neigbor=1,2

          ENDIF !  patch_3D%lsm_e(je,jk,block) <= sea_boundary
          
        END DO
      END DO
    END DO

    DO jev=1,2*no_dual_edges
      DO jk = 1, n_zlev
        CALL sync_patch_array(SYNC_E, patch_2D, operators_coefficients%edge2edge_viavert_coeff(:,jk,:, jev))
      ENDDO
    ENDDO

    !-------------------------------------------------------------

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE apply_boundary2coeffs
!-------------------------------------------------------------------------
  !>
  !! Precomputes the geometrical factors used in the divergence, rotation.
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE init_diff_operator_coeff_3D( patch_2D, operators_coefficients )
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: patch_2D
    TYPE(t_operator_coeff),     INTENT(inout) :: operators_coefficients
    !

    INTEGER :: ie
    !INTEGER :: rl_start, rl_end
    INTEGER :: i_nchdom!,i_startblk, i_endblk, i_startidx, i_endidx
  
    TYPE(t_subset_range), POINTER :: all_edges, owned_cells, owned_verts
    INTEGER :: edge_block, edge_index
    INTEGER :: cell_index, cell_block
    INTEGER :: vertex_index, vertex_block
    INTEGER :: start_index, end_index, neigbor, level

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_operator_ocean_coeff_3d:init_diff_operator_coeff_3D')
    !-----------------------------------------------------------------------
    i_nchdom   = MAX(1,patch_2D%n_childdom)
    all_edges => patch_2D%edges%all
    owned_cells => patch_2D%cells%owned
    owned_verts => patch_2D%verts%owned

!!$OMP PARALLEL
    ! 1) coefficients for divergence
!!$OMP DO PRIVATE(cell_block, start_index, end_index, cell_index, neigbor, &
!!$OMP edge_index, edge_block, level) ICON_OMP_DEFAULT_SCHEDULE
    DO cell_block = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, cell_block, start_index, end_index)
      DO cell_index = start_index, end_index
        DO neigbor=1, patch_2D%cells%num_edges(cell_index, cell_block)

          edge_index = patch_2D%cells%edge_idx(cell_index, cell_block, neigbor)
          edge_block = patch_2D%cells%edge_blk(cell_index, cell_block, neigbor)

          operators_coefficients%div_coeff(cell_index, 1, cell_block, neigbor) = &
            & patch_2D%edges%primal_edge_length(edge_index, edge_block)        * &
            & patch_2D%cells%edge_orientation(cell_index, cell_block, neigbor) / &
            & patch_2D%cells%area(cell_index, cell_block)

          DO level=2, n_zlev
            operators_coefficients%div_coeff(cell_index, level, cell_block, neigbor) = &
               operators_coefficients%div_coeff(cell_index, 1, cell_block, neigbor)
          ENDDO !levels

        ENDDO !neigbor
      ENDDO ! cell_index = start_index, end_index
    ENDDO !cell_block
!!$OMP ENDDO NOWAIT

    ! 2) coefficients for curl
!!$OMP DO PRIVATE(vertex_block, start_index, end_index, vertex_index, neigbor, &
!!$OMP edge_index, edge_block, level) ICON_OMP_DEFAULT_SCHEDULE
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index
        DO neigbor=1, patch_2D%verts%num_edges(vertex_index, vertex_block)
          edge_index = patch_2D%verts%edge_idx(vertex_index, vertex_block, neigbor)
          edge_block = patch_2D%verts%edge_blk(vertex_index, vertex_block, neigbor)

          operators_coefficients%rot_coeff(vertex_index,1,vertex_block,neigbor)  =     &
            & dual_edge_length(edge_index,edge_block) * grid_sphere_radius *     &
            & patch_2D%verts%edge_orientation(vertex_index,vertex_block,neigbor)

          DO level=2, n_zlev
            operators_coefficients%rot_coeff(vertex_index,level,vertex_block,neigbor) = &
               operators_coefficients%rot_coeff(vertex_index,1,vertex_block,neigbor)
          ENDDO !levels
        
        ENDDO !neigbor
      ENDDO ! vertex_index = start_index, end_index
    ENDDO !vertex_block = owned_verts%start_block, owned_verts%end_block
!!$OMP END DO NOWAIT


    ! 4) coefficients for gradient
!!$OMP DO PRIVATE(edge_block, start_index, end_index, edge_index, level) ICON_OMP_DEFAULT_SCHEDULE
    DO edge_block = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, edge_block, start_index, end_index)
      DO edge_index = start_index, end_index
        ! this should be calculated in the patch_2D setup
        patch_2D%edges%inv_dual_edge_length(edge_index,edge_block) = &
          & 1._wp / (dual_edge_length(edge_index,edge_block)  * grid_sphere_radius)

        DO level=1, n_zlev
          operators_coefficients%grad_coeff(edge_index,level, edge_block)   &
            & = patch_2D%edges%inv_dual_edge_length(edge_index,edge_block)
        ENDDO !levels

      ENDDO ! edge_index = start_index, end_index
    ENDDO  ! edge_block
!!$OMP ENDDO NOWAIT
!!$OMP END PARALLEL

    ! no need to synchronize all elements of operators_coefficients%grad_coeff
!     CALL sync_patch_array(sync_e, patch_2D, operators_coefficients%grad_coeff)

    DO ie = 1, no_primal_edges
      CALL sync_patch_array(sync_c, patch_2D, operators_coefficients%div_coeff(:,:,:,ie))
    END DO

    DO ie = 1, no_dual_edges
      CALL sync_patch_array(sync_v, patch_2D, operators_coefficients%rot_coeff(:,:,:,ie))
    END DO

  END SUBROUTINE init_diff_operator_coeff_3D
  !-------------------------------------------------------------------------


END MODULE mo_operator_ocean_coeff_3d

