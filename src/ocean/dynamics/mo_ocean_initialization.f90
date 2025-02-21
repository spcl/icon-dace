!  Contains the data structures for the hydrostatic ocean model.
!
!  Contains the data structures to store the hydrostatic & boussinesq ocean model state.
!  Implementation is based on ICON-Shallow-Water model
!  to store the shallow water model state and other auxiliary variables.
!  Constructors and destructors for these data structures are also defined here.
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
MODULE mo_ocean_initialization
  !-------------------------------------------------------------------------
  USE mo_kind,                ONLY: wp
  USE mo_mpi,                 ONLY: my_process_is_mpi_test, p_pe
  USE mo_parallel_config,     ONLY: nproma
  USE mo_master_config,       ONLY: isRestart
  USE mo_impl_constants,      ONLY: land, land_boundary, boundary, sea_boundary, sea,  &
    & success, max_char_length, min_dolic,               &
    & full_coriolis, beta_plane_coriolis,                &
    & f_plane_coriolis, zero_coriolis, halo_levels_ceiling, &
    & on_cells, on_edges, on_vertices
  USE mo_ocean_nml,           ONLY: n_zlev, dzlev_m, no_tracer, l_max_bottom, l_partial_cells, &
    & coriolis_type, basin_center_lat, basin_height_deg, iswm_oce, coriolis_fplane_latitude,   &
    & use_smooth_ocean_boundary, max_allocated_levels, minVerticalLevels
  USE mo_util_dbg_prnt,       ONLY: c_i, c_b, nc_i, nc_b, near_proc_id
  USE mo_exception,           ONLY: message_text, message, finish
  USE mo_model_domain,        ONLY: t_patch,t_patch_3d, t_grid_cells, t_grid_edges
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, grid_sphere_radius, grid_angular_velocity, &
    & use_dummy_cell_closure
  USE mo_ext_data_types,      ONLY: t_external_data
  USE mo_dynamics_config,     ONLY: nnew,nold
  USE mo_math_types,          ONLY: t_cartesian_coordinates, t_geographical_coordinates
  USE mo_math_utilities,      ONLY: gc2cc, arc_length, set_zlev
  USE mo_math_constants,      ONLY: deg2rad,rad2deg
  USE mo_sync,                ONLY: sync_e, sync_c, sync_v,sync_patch_array, global_sum_array, sync_idx, &
    & enable_sync_checks, disable_sync_checks
  USE mo_cf_convention
!  USE mo_grib2
  USE mo_grid_subset,         ONLY: t_subset_range, get_index_range, fill_subset
  ! USE mo_ocean_config,        ONLY: ignore_land_points
  USE mo_ocean_types, ONLY: t_hydro_ocean_base
  USE mo_ocean_diagnostics_types, ONLY: &
    & t_ocean_regions, &
    & t_ocean_region_volumes, &
    & t_ocean_region_areas, &
    & t_ocean_basins
  USE mo_ocean_state, ONLY:  ocean_restart_list, &
    & ocean_default_list, &
    & v_base
  USE mo_util_dbg_prnt,       ONLY: dbg_print, debug_print_MaxMinMean

  USE mo_ocean_check_tools, ONLY: ocean_check_level_sea_land_mask, check_ocean_subsets

  USE mo_mpi,                ONLY: my_process_is_stdio, get_my_mpi_work_id

  IMPLICIT NONE
  PRIVATE


  !public interface
  !
  ! subroutines
  ! PUBLIC :: set_lateral_boundary_values
  PUBLIC :: init_ho_base
  PUBLIC :: init_ho_basins
  PUBLIC :: init_coriolis_oce
  PUBLIC :: is_initial_timestep
  PUBLIC :: check_ocean_subsets

  PUBLIC :: init_patch_3d
  !

CONTAINS




  !-------------------------------------------------------------------------
  !>
  !! Sbr set boundary values for velocity field.
  !!
  !
!  SUBROUTINE set_lateral_boundary_values( patch_3d, vn )
!
!    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
!    REAL(wp)                             :: vn(:,:,:)
!
!    ! local variables
!    INTEGER :: jb, je, jk
!    INTEGER :: StartEdgeIndex, EndEdgeIndex
!    INTEGER :: slev,elev
!    TYPE(t_subset_range), POINTER :: all_edges
!    TYPE(t_patch), POINTER :: patch_2d
!    !---------------------------------------------------------------
!    patch_2d   => patch_3d%p_patch_2d(1)
!    all_edges => patch_2d%edges%ALL
!
!    ! blocking
!    slev         = 1
!    elev         = n_zlev
!
!    DO jb = all_edges%start_block, all_edges%end_block
!      CALL get_index_range(all_edges, jb, StartEdgeIndex, EndEdgeIndex)
!      DO jk = slev, elev
!        DO je= StartEdgeIndex, EndEdgeIndex
!          IF ( patch_3d%lsm_e(je,jk,jb) >= boundary ) THEN
!            vn(je,jk,jb) = 0.0_wp
!          ENDIF
!        END DO
!      END DO
!    END DO
!  END SUBROUTINE set_lateral_boundary_values
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Initializes 3-dimensional structures for the ocean state variable.
  !!
  !! Calculates the vertical grid levels in z (meter) using the thickness of the elemental
  !! prisms (del_zlev_i) that are read from the ocean namelist.
  !!
  !! The 3-dimensional land-sea-mask is filled with values for interieur ocean, boundary
  !! ocean, and land, where parameter values from mo_impl_constants are used. The three
  !! dimensions are two for the nproma-blocking  and the middle one for the vertical levels.
  !!
  !!
  !!  mpi parallelized
!<Optimize:inUse>
  SUBROUTINE init_ho_base( patch_2d, p_ext_data, v_base )

    TYPE(t_patch),  TARGET,   INTENT(inout)    :: patch_2d
    TYPE(t_external_data),    INTENT(inout)    :: p_ext_data
    TYPE(t_hydro_ocean_base), INTENT(inout)    :: v_base

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_initialization:init_ho_base'

    INTEGER :: jb, jc, je, jk, ji, alloc_cell_blocks, nblks_e, npromz_c, npromz_e
    INTEGER :: i_startidx, i_endidx
    INTEGER :: noct1_c, noct1_e, inolsm, nowet_c
    INTEGER :: nolnd_c(n_zlev), nosea_c(n_zlev), nogllnd_c, noglsea_c
    INTEGER :: nolnd_e(n_zlev), nosea_e(n_zlev), nogllnd_e, noglsea_e
    INTEGER :: nobnd_e(n_zlev), nosbd_c(n_zlev), nolbd_c(n_zlev)
    INTEGER :: noglbnd_e, noglsbd_c, nogllbd_c
    INTEGER :: iic1, ibc1, iic2, ibc2, idxe, ible, idxn, ibln
    INTEGER :: n_zlvp, n_zlvm
    INTEGER :: jiter, niter, ctr, ctr_jk, ctr_glb
    INTEGER :: lsm_c   (nproma,patch_2d%alloc_cell_blocks)
    REAL(wp):: perc_lnd_c(n_zlev), perc_gllnd_c
    REAL(wp):: perc_lnd_e(n_zlev), perc_gllnd_e
    REAL(wp):: z_stepvalue

    REAL(wp):: z_sync_c(nproma,patch_2d%alloc_cell_blocks)
    REAL(wp):: z_sync_e(nproma,patch_2d%nblks_e)
    REAL(wp):: z_lat, z_lat_deg, z_north, z_south

    TYPE(t_subset_range), POINTER :: owned_cells, all_cells
    TYPE(t_subset_range), POINTER :: edges_in_domain, all_edges
    INTEGER :: all_nobnd_e, all_nosbd_c, all_nolbd_c
    INTEGER :: dol_e, dol_c1, dol_c2, lsm_e

    LOGICAL :: limited_area
    LOGICAL :: l_vert_step

    limited_area = .FALSE.
    l_vert_step  = .FALSE.
    !-----------------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')

    owned_cells => patch_2d%cells%owned
    all_cells   => patch_2d%cells%ALL
    edges_in_domain => patch_2d%edges%in_domain
    all_edges   => patch_2d%edges%ALL

    z_sync_c(:,:) = 0.0_wp
    z_sync_e(:,:) = 0.0_wp
    lsm_c   (:,:) = 0

    z_lat     = 0.0_wp
    z_lat_deg = 0.0_wp
    z_north   = 0.0_wp
    z_south   = 0.0_wp

    !-----------------------------
    !
    ! Basic z-level configuration:
    !
    ! n_zlev    : number of z-coordinate surfaces - module global_variables, read in from namelist
    ! del_zlev_m: thickness of elemental prism - read in from namelist
    ! zlev_m    : position of coordinate surfaces in meters below zero surface
    ! zlev_i    : surface at the top of the respective z-coordinate surface (intermediate level)
    ! del_zlev_i: distance between two z-coordinate surfaces
    !
    !-----------------------------

    ! values for the blocking
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    nblks_e = patch_2d%nblks_e
    !nblks_v = patch_2D%nblks_v
    npromz_c = patch_2d%npromz_c
    npromz_e = patch_2d%npromz_e

    ! number of vertical levels from the namelist in mo_global_variables
    n_zlvp = n_zlev + 1
    n_zlvm = n_zlev - 1

    !-----------------------------
    ! Create a step inside the deep ocean bathymetry
    !-----------------------------

    IF (l_vert_step .AND. p_pe == near_proc_id) THEN

      ! Set the bathymetry at cell and edges to a value for tests mainly with Aqua Planet

      !  - blocks and indices of dbg_print are used, i.e. location of the step
      !    can be controled by dbg_lat/lon_in of namelist dbg_index_nml
      !  - Attention: lsm for level 1 and 2 is set by gridgenerator, therefore minimum depth
      !    must be greater than zlev_i(3)
      !  - Attention: when running in parallel there might be steps as much as mpi-processes
      !    or even errors - use in serial mode only
      !  - Attention: bathymetry is negative in the ocean

      z_stepvalue = -100.0_wp
      p_ext_data%oce%bathymetry_c(c_i,c_b) = z_stepvalue
      !  first neighbor gets another value:
      p_ext_data%oce%bathymetry_c(nc_i(1),nc_b(1)) = z_stepvalue - 520.0_wp

    END IF

    !-----------------------------
    !
    ! Fill the 3-dim land-sea-mask and number of deepest ocean layer in column 'dolic'
    !
    !-----------------------------

    ! fill one-dimensional vertical levels
    CALL set_del_zlev(n_zlev, dzlev_m,                  &
      & v_base%del_zlev_i, v_base%del_zlev_m, &
      & v_base%zlev_i    , v_base%zlev_m)


    nogllnd_c = 0
    noglsea_c = 0
    nogllnd_e = 0
    noglsea_e = 0
    noglbnd_e = 0
    noglsbd_c = 0
    nogllbd_c = 0
    noct1_c = 0
    noct1_e = 0

    ! surface level: as read in ext_data
    ! for cells fill all levels in orddr to take care of ghost land cells
    DO jk = 1, n_zlev
      v_base%lsm_c(:,jk,:) = p_ext_data%oce%lsm_ctr_c(:,:)
    ENDDO
    v_base%lsm_e(:,1,:) = p_ext_data%oce%lsm_ctr_e(:,:)
    !     !  surface level and second level of lsm_c defined by gridgenerator, not the current bathymetry
    !     v_base%lsm_c(:,1,:) = p_ext_data%oce%lsm_ctr_c(:,:)
    !     IF(n_zlev>=2) v_base%lsm_c(:,2,:) = p_ext_data%oce%lsm_ctr_c(:,:)

    !  first and second level of dolic_c defined by gridgenerator
    WHERE (p_ext_data%oce%lsm_ctr_c(:,:) <= sea_boundary) v_base%dolic_c(:,:) = 2


    ! Coordinate surfaces - n_zlev z-levels:
    ! First vertical level loop to set wet cells below surface (and second layer) only

    init_zloop: DO jk = minVerticalLevels+1, n_zlev

      !-----------------------------
      ! set dolic and wet grid points on cells:
      !  - if bathymetry is deeper than or equal to the top of the coordinate surface,
      !    (zlev_i(jk)), then grid point is wet; dolic is in that level (l_max_bottom=.true.)
      !  - if bathymetry is deeper than or equal to the coordinate surface (zlev_m)
      !    then grid point is wet; dolic is in that level (l_max_bottom=.false.)
      !  - values for BOUNDARY set below

!      IF (l_partial_cells) l_max_bottom = .FALSE.
      IF (l_max_bottom) THEN

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
          DO jc = i_startidx, i_endidx

            IF (p_ext_data%oce%bathymetry_c(jc,jb) <= -v_base%zlev_i(jk)) THEN
              v_base%lsm_c(jc,jk,jb) = sea
              v_base%dolic_c(jc,jb)  = jk
            ELSE
              v_base%lsm_c(jc,jk,jb) = land
            END IF

          END DO
        END DO

      ELSE

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(all_cells, jb, i_startidx, i_endidx)
          DO jc = i_startidx, i_endidx

            IF (p_ext_data%oce%bathymetry_c(jc,jb) <= -v_base%zlev_m(jk)) THEN
              v_base%lsm_c(jc,jk,jb) = sea
              v_base%dolic_c(jc,jb)  = jk
            ELSE
              v_base%lsm_c(jc,jk,jb) = land
            END IF

          END DO
        END DO

      END IF

      ! synchronize lsm on cells
      ! LL: this is done on all cells consistently on al procs
      ! so no sync is required

    END DO init_zloop


    !-------------------------------------------------
    IF(limited_area)THEN
      !  #slo# must be parallelized accordingly

      DO jk = 1, n_zlev
        z_south=-80.0_wp

        DO jb = all_cells%start_block, all_cells%end_block
          CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)
          DO jc = i_startidx, i_endidx

            !get latitude of actual cell
            z_lat = patch_2d%cells%center(jc,jb)%lat
            z_lat_deg = z_lat*rad2deg

            !If latitude of cell is above 80 N or below 80 S set triangle to land
            IF(z_lat_deg>z_north.OR.z_lat_deg<z_south)THEN
              v_base%lsm_c(jc,:,jb)          = land
              p_ext_data%oce%bathymetry_c(jc,jb) = 100.0_wp
              v_base%dolic_c(jc,jb)              = 0
              v_base%wet_c(jc,:,jb)              = 0.0_wp
              !Set also all 3 edges to land
              DO ji = 1, 3
                ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
                idxe                                   = patch_2d%cells%edge_idx(jc,jb,ji)
                ible                                   = patch_2d%cells%edge_blk(jc,jb,ji)
                v_base%lsm_e(idxe,:,ible)          = land
                p_ext_data%oce%bathymetry_e(idxe,ible) = 100.0_wp
                v_base%dolic_e(idxe,ible)              = 0
                v_base%wet_e(idxe,:,ible)              = 0.0_wp
              END DO
            ENDIF
          END DO
        END DO
      END DO
    ENDIF

    !---------------------------------------------------------------------------------------------
    ! Correction loop for cells in all levels, similar to surface done in grid generator
    !  - through all levels each wet cell has at most one dry cell as neighbour
    !  - otherwise it is set to dry grid cell
    IF (use_smooth_ocean_boundary) THEN
      niter=30
      zloop_cor: DO jk=1,n_zlev

        ctr_jk = 0

        ! working on 2D lsm_c inside the loop
        lsm_c(:,:) = v_base%lsm_c(:,jk,:)

        ! LL: disable checks here, the changes in halos will differ from seq run
        !     as the access patterns differ
        CALL disable_sync_checks()

        DO jiter=1,niter
          !
          ctr = 0 ! no changes initially

          ! loop through owned patch cells
          DO jb = owned_cells%start_block, owned_cells%end_block
            CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)

            DO jc =  i_startidx, i_endidx

              nowet_c = 0

              ! LL: here we probably want to check if the above cell is land
              !     and change this into land accordingly

              IF (lsm_c(jc,jb) <= sea_boundary) THEN
                DO ji = 1, 3
                  ! Get indices/blks of cells 1 to 3 adjacent to cell (jc,jb)
                  idxn = patch_2d%cells%neighbor_idx(jc,jb,ji)
                  ibln = patch_2d%cells%neighbor_blk(jc,jb,ji)
                  ! counts number of land-cells for all three neighbors
                  !  - only one land-point neighbor is allowed
                  IF (idxn <= 0) THEN
                    nowet_c = nowet_c + 1
                  ELSE IF ( lsm_c(idxn,ibln) > sea_boundary ) THEN
                    nowet_c = nowet_c + 1
                  ENDIF
                END DO

                ! More than 1 wet neighbor-cell then set cell to land
                !  - edges are set in the correction loop below
                IF ( nowet_c >= 2 ) THEN
                  lsm_c(jc,jb)=land_boundary
                  ctr = ctr+1

                  IF (jk<3) THEN
                    WRITE(message_text,'(a,2i8)') &
                      & 'WARNING: Found 2 land neighbors at jc, jk=',jc,jk
                    CALL message(TRIM(routine), TRIM(message_text))
                  END IF

                END IF ! 2 land neighbors

              END IF ! lsm_c(jc,jb) <= SEA_BOUNDARY

            END DO  ! jc =  i_startidx, i_endidx
          END DO ! jb = owned_cells%start_block, owned_cells%end_block

          ! see what is the sum of changes of all procs
          ctr_glb = global_sum_array(ctr)

          WRITE(message_text,'(a,i2,a,i2,a,i8)') 'Level:', jk, &
            & ' Corrected wet cells with 2 land neighbors - iter=', &
            & jiter,' no of cor:',ctr_glb
          CALL message(TRIM(routine), TRIM(message_text))

          ! if no changes have been done, we are done with this level. Exit
          IF (ctr_glb == 0) EXIT

          ! we need to sync the halos here
          z_sync_c(:,:) =  REAL(lsm_c(:,:),wp)
          CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
          lsm_c(:,:) = INT(z_sync_c(:,:))

        END DO   ! jiter

        CALL enable_sync_checks()

        z_sync_c(:,:) = REAL(lsm_c(:,:),wp)
        CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
        lsm_c(:,:) = INT(z_sync_c(:,:))

        ! get back into 3D the slm
        v_base%lsm_c(:,jk,:) = lsm_c(:,:)

      END DO zloop_cor ! jk=1,n_zlev
    ENDIF !(use_smooth_ocean_boundary)

    ! restore p_test_run
    CALL enable_sync_checks()

    !---------------------------------------------------------------------------------------------
    ! Now run through whole zlevel_loop after correction of cells for calculation of boundaries

    nogllnd_c = 0
    noglsea_c = 0
    nogllnd_e = 0
    noglsea_e = 0
    noglbnd_e = 0
    noglsbd_c = 0
    nogllbd_c = 0

    ! set dolic once more after jiter-correction
    v_base%dolic_c = 0
    v_base%dolic_e = 0

    ! Main loop for edges, dolic, boundaries, diagnosis and output
    !  - using lsm_c after complete correction in jk>2 as input
    !  - (1) set land and sea values at cells to SEA=-2 and LAND=2 without considering boundaries
    !  - (2) set land and sea values at edges including boundaries
    !  - (3) set land and sea boundary values (-1 = SEA_BOUNDARY, 1=LAND_BOUNDARY)
    !
    zlevel_loop: DO jk = 1, n_zlev

      !-----------------------------
      ! (1) set wet grid points and dolic at cells:
      !  - values for BOUNDARY set below

      nolnd_c(jk)=0
      nosea_c(jk)=0


      DO jb = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)
        DO jc = i_startidx, i_endidx

          ! IF (.NOT.patch_2D%cells%decomp_info%owner_mask(jc,jb)) CYCLE  ! access inner domain only
          IF (v_base%lsm_c(jc,jk,jb) <= sea_boundary) THEN
            nosea_c(jk)=nosea_c(jk)+1
            v_base%dolic_c(jc,jb) = jk
          ELSE
            ! -after correction: all other grid points are set to dry
            v_base%lsm_c(jc,jk,jb) = land
            nolnd_c(jk)=nolnd_c(jk)+1
          END IF

          ! counting surface conditions as read from bathymetry - all wet cells
          IF (jk == 1) THEN
            IF (p_ext_data%oce%bathymetry_c(jc,jb) <= -v_base%zlev_m(jk)) &
              & noct1_c = noct1_c+1
          ENDIF

        END DO

      END DO

      ! now synchronize lsm_c
      z_sync_c(:,:) =  REAL(v_base%lsm_c(:,jk,:),wp)
      CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
      v_base%lsm_c(:,jk,:) = INT(z_sync_c(:,:))

      !  percentage of land area per level and global value
      !   - here: nosea/nolnd include boundaries
      ctr         = global_sum_array(nolnd_c(jk))
      nolnd_c(jk) = ctr
      ctr         = global_sum_array(nosea_c(jk))
      nosea_c(jk) = ctr
      inolsm      = nolnd_c(jk) + nosea_c(jk)
      IF (inolsm == 0 ) THEN
        IF (jk < 3 ) CALL message (TRIM(routine), 'WARNING - number of cell points is zero?')
        perc_lnd_c(jk) = 0.0_wp
      ELSE
        perc_lnd_c(jk) = REAL(nolnd_c(jk),wp)/REAL(nosea_c(jk)+nolnd_c(jk),wp)*100.0_wp
        nogllnd_c = nogllnd_c + nolnd_c(jk)
        noglsea_c = noglsea_c + nosea_c(jk)
      END IF


      !-----------------------------
      ! (2) set wet grid points and dolic at edges (get values of neighbouring cells)
      !  - if the two corresponding cells are differing in sign then edge is BOUNDARY
      !  - if the two corresponding cells are <0 then edge is SEA
      !  - if the two corresponding cells are >0 then edge is LAND

      nolnd_e(jk)=0
      nosea_e(jk)=0
      nobnd_e(jk)=0

      ! loop through owned patch edges
      DO jb = edges_in_domain%start_block, edges_in_domain%end_block
        CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
        DO je = i_startidx, i_endidx

          ! IF (.NOT. patch_2D%edges%decomp_info%owner_mask(je,jb)) CYCLE  ! access inner domain only

          ! get indices/blks of cells 1 and 2 adjacent to edge (je,jb)
          iic1 = patch_2d%edges%cell_idx(je,jb,1)
          ibc1 = patch_2d%edges%cell_blk(je,jb,1)
          iic2 = patch_2d%edges%cell_idx(je,jb,2)
          ibc2 = patch_2d%edges%cell_blk(je,jb,2)

          ! when a cell is missing then the edge is considered boundary (=0)
          ! if the other cell is sea
          ! in case of open boundaries we have to create a better process through
          ! the grid tools
          IF (iic1 == 0 ) THEN
            IF (v_base%lsm_c(iic2,jk,ibc2) <= sea_boundary) THEN
              v_base%lsm_e(je,jk,jb) = boundary
            ELSE
              v_base%lsm_e(je,jk,jb) = land
            ENDIF
          ELSE IF (iic2 == 0 ) THEN
            IF (v_base%lsm_c(iic1,jk,ibc1) <= sea_boundary) THEN
              v_base%lsm_e(je,jk,jb) = boundary
            ELSE
              v_base%lsm_e(je,jk,jb) = land
            ENDIF
          ELSE

            ! set land/sea for all edges
            IF ( (v_base%lsm_c(iic1,jk,ibc1) <= sea_boundary)  .AND.   &
              & (v_base%lsm_c(iic2,jk,ibc2) <= sea_boundary) ) THEN
              v_base%lsm_e(je,jk,jb) = sea
            ELSEIF ( (v_base%lsm_c(iic1,jk,ibc1) >= LAND_BOUNDARY)  .AND.   &
              & (v_base%lsm_c(iic2,jk,ibc2) >= LAND_BOUNDARY) ) THEN
              v_base%lsm_e(je,jk,jb) = land
            ELSE
              v_base%lsm_e(je,jk,jb) = boundary
            ENDIF

          ENDIF

          IF ( patch_2d%edges%decomp_info%owner_mask(je,jb)) THEN
            ! count land/sea/boundary values (sum of nosea_e no_lnd_e nobnd_e is global value)
            IF ( v_base%lsm_e(je,jk,jb) <  boundary )      &
              & nosea_e(jk)=nosea_e(jk)+1
            IF ( v_base%lsm_e(je,jk,jb) >  boundary )      &
              & nolnd_e(jk)=nolnd_e(jk)+1
            IF ( v_base%lsm_e(je,jk,jb) == boundary )      &
              & nobnd_e(jk)=nobnd_e(jk)+1
          ENDIF

          ! set dolic to jk if lsm_e is wet (minimum depth of 2 neighboring cells)
          IF ( v_base%lsm_e(je,jk,jb) < boundary )      &
            & v_base%dolic_e(je,jb) = jk

          ! counting surface conditions as read from bathymetry - all wet edges, boundary edges
          IF (jk == 1) THEN
            IF (p_ext_data%oce%bathymetry_e(je,jb) <= -v_base%zlev_m(jk)) noct1_e = noct1_e+1
          ENDIF

        END DO

      END DO

      ! synchronize lsm on edges
      z_sync_e(:,:) =  REAL(v_base%lsm_e(:,jk,:),wp)
      CALL sync_patch_array(sync_e, patch_2d, z_sync_e(:,:))
      v_base%lsm_e(:,jk,:) = INT(z_sync_e(:,:))

      !  percentage of land area per level and global value
      ctr         = global_sum_array(nolnd_e(jk))
      nolnd_e(jk) = ctr
      ctr         = global_sum_array(nosea_e(jk))
      nosea_e(jk) = ctr
      inolsm = nolnd_e(jk) + nosea_e(jk)
      IF (inolsm == 0 ) THEN
        IF (jk == 1 ) CALL message (TRIM(routine), 'WARNING - number of edge points is zero?')
        perc_lnd_e(jk) = 0.0_wp
      ELSE
        perc_lnd_e(jk) = REAL(nolnd_e(jk),wp)/REAL(nosea_e(jk)+nolnd_e(jk),wp)*100.0_wp
        nogllnd_e = nogllnd_e + nolnd_e(jk)
        noglsea_e = noglsea_e + nosea_e(jk)
      END IF

      !-----------------------------
      ! (3) set values for LAND_BOUNDARY and SEA_BOUNDARY at cells
      !  - get values of neighbouring edges
      !  - if 1 of 3 edges of a sea-cell is BOUNDARY then cell is SEA_BOUNDARY
      !  - if 1 (or 2) of 3 edges of a land-cell is BOUNDARY then cell is LAND_BOUNDARY

      nosbd_c(jk)=0
      nolbd_c(jk)=0

      ! loop through owned patch cells since we count changes globally
      DO jb = owned_cells%start_block, owned_cells%end_block
        CALL get_index_range(owned_cells, jb, i_startidx, i_endidx)
        DO jc =  i_startidx, i_endidx

          ! sea points
          IF (v_base%lsm_c(jc,jk,jb) < boundary) THEN

            DO ji = 1, 3
              ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
              idxe = patch_2d%cells%edge_idx(jc,jb,ji)
              ible = patch_2d%cells%edge_blk(jc,jb,ji)
              ! if one of lsm_e is boundary then lsm_c is sea_boundary
              IF ( v_base%lsm_e(idxe,jk,ible) == boundary ) &
                & v_base%lsm_c(jc,jk,jb) = sea_boundary
            END DO

            ! count sea boundary for all levels
            IF ( v_base%lsm_c(jc,jk,jb) == sea_boundary )  &
              & nosbd_c(jk)=nosbd_c(jk)+1
          END IF  !  lsm_c < 0

          ! land points
          IF (v_base%lsm_c(jc,jk,jb) > boundary) THEN

            DO ji = 1, 3
              ! Get indices/blks of edges 1 to 3 adjacent to cell (jc,jb)
              idxe = patch_2d%cells%edge_idx(jc,jb,ji)
              ible = patch_2d%cells%edge_blk(jc,jb,ji)
              ! if one of lsm_e is boundary then lsm_c is land_boundary
              IF ( v_base%lsm_e(idxe,jk,ible) == boundary ) &
                & v_base%lsm_c(jc,jk,jb) = land_boundary
            END DO

            IF ( v_base%lsm_c(jc,jk,jb) == land_boundary )   &
              & nolbd_c(jk)=nolbd_c(jk)+1

          END IF  !  lsm_c > 0

        END DO

      END DO

      all_nobnd_e = global_sum_array( nobnd_e(jk))
      all_nosbd_c = global_sum_array( nosbd_c(jk))
      all_nolbd_c = global_sum_array( nolbd_c(jk))

      nobnd_e(jk) = all_nobnd_e
      nosbd_c(jk) = all_nosbd_c
      nolbd_c(jk) = all_nolbd_c

      noglbnd_e = noglbnd_e + all_nobnd_e
      noglsbd_c = noglsbd_c + all_nolbd_c
      nogllbd_c = nogllbd_c + all_nolbd_c

      ! synchronize lsm on cells
      z_sync_c(:,:) =  REAL(v_base%lsm_c(:,jk,:),wp)
      CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
      v_base%lsm_c(:,jk,:) = INT(z_sync_c(:,:))

    END DO zlevel_loop

    ctr     = global_sum_array( noct1_c )
    noct1_c = ctr
    ctr     = global_sum_array( noct1_e )
    noct1_e = ctr

    !---------------------------------------------------------------------------------------------
    ! synchronize dolic_c
    z_sync_c(:,:) =  REAL(v_base%dolic_c(:,:),wp)
    CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
    v_base%dolic_c(:,:) = INT(z_sync_c(:,:))
    ! synchronize dolic_e
    z_sync_e(:,:) =  REAL(v_base%dolic_e(:,:),wp)
    CALL sync_patch_array(sync_e, patch_2d, z_sync_e(:,:))
    v_base%dolic_e(:,:) = INT(z_sync_e(:,:))

    !---------------------------------------------------------------------------------------------
    ! Output the levels
    WRITE(message_text,'(a,a)') &
      & 'LEVEL   zlev_m  Thickness   zlev_i  Distance ', &
      & '    SEA_c    LAND_c  PERC_LND     SEA_e    LAND_e   BND_e   SEA_B. LAND_B.'
    CALL message('', TRIM(message_text))
    DO jk = 1, n_zlev
      WRITE(message_text,'(a,i3,4f10.2,2i10,f10.2,2i10,3i8)') '.',  &
        & jk, v_base%zlev_m(jk), v_base%del_zlev_m(jk), &
        & v_base%zlev_i(jk), v_base%del_zlev_i(jk), &
        & nosea_c(jk), nolnd_c(jk), perc_lnd_c(jk), nosea_e(jk), nolnd_e(jk),     &
        & nobnd_e(jk), nosbd_c(jk), nolbd_c(jk)
      CALL message('', message_text)
    END DO

    ! Output last level
    inolsm = nogllnd_c + noglsea_c
    IF ( inolsm == 0 ) THEN
      CALL message (TRIM(routine), 'WARNING - number of global cell points is zero?')
      perc_gllnd_c = 0.0_wp
    ELSE
      perc_gllnd_c = REAL(nogllnd_c,wp)/REAL(noglsea_c+nogllnd_c,wp)*100.0_wp
    END IF
    inolsm = nogllnd_e + noglsea_e
    IF ( inolsm == 0 ) THEN
      CALL message (TRIM(routine), 'WARNING - number of global edge points is zero?')
      perc_gllnd_e = 0.0_wp
    ELSE
      perc_gllnd_e = REAL(nogllnd_e,wp)/REAL(noglsea_e+nogllnd_e,wp)*100.0_wp
    END IF
    n_zlvp = n_zlev + 1
    WRITE(message_text,'(a,f20.2,a,i9,i10,f10.2,2i10,3i8)') &
      & 'Bottom Level: ', v_base%zlev_i(n_zlvp), &
      & '    GLOBAL:',     noglsea_c, nogllnd_c, perc_gllnd_c, &
      & noglsea_e, nogllnd_e, noglbnd_e, noglsbd_c, nogllbd_c
    CALL message('', TRIM(message_text))

    !---------------------------------------------------------------------------------------------
    ! CHECKS

    ! Warnings occur if create_ocean_grid parameter mindepth is not half the
    ! depth of the surface level dzlev_m(1) in namelist ocean_ctl (must not be an error)
    IF ( nosea_c(1) /= noct1_c ) THEN
      WRITE(message_text,'(a,i8,a,i8)') &
        & 'WARNING - surface sea-cells read = ',nosea_c(1), &
        & ' - calculated from bathymetry = ',noct1_c
      CALL message(routine, TRIM(message_text))
    END IF
    IF ( nosea_e(1) /= noct1_e ) THEN
      WRITE(message_text,'(a,i8,a,i8)') &
        & 'WARNING - surface sea-edges read = ',nosea_e(1), &
        & ' - calculated from bathymetry = ',noct1_e
      CALL message(routine, TRIM(message_text))
    END IF

    ! Test loop for dolic_e:
    !  - if lsm_e(dolic_e+1) is boundary, then dolic_c1 or c2 > dolic_e
    !  - more tests? lsm(dolic) is no boundary any more
    !  - bugfix: edges_in_domain for test only
    !TODO: review usage of v_base, owned vs. all_edges
    ! this is done in ocean_check_level_sea_land_mask
    !    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
    !      CALL get_index_range(edges_in_domain, jb, i_startidx, i_endidx)
    !      DO je = i_startidx, i_endidx
    !
    !        dol_e = v_base%dolic_e(je,jb)
    !
    !        IF (dol_e > 1 .AND. dol_e < n_zlev) THEN
    !
    !          lsm_e = v_base%lsm_e(je,dol_e+1,jb)
    !
    !          ! get indices/blks of cells 1 and 2 adjacent to edge (je,jb)
    !          iic1 = patch_2D%edges%cell_idx(je,jb,1)
    !          ibc1 = patch_2D%edges%cell_blk(je,jb,1)
    !          iic2 = patch_2D%edges%cell_idx(je,jb,2)
    !          ibc2 = patch_2D%edges%cell_blk(je,jb,2)
    !
    !          IF (iic1 > 0 .AND. iic2 > 0) THEN
    !            dol_c1 = v_base%dolic_c(iic1,ibc1)
    !            dol_c2 = v_base%dolic_c(iic2,ibc2)
    !
    !            IF (dol_c1 == dol_c2 .AND. lsm_e == 0) THEN
    !              WRITE(message_text,'(a,2i3,a,i3)') &
    !                &   'WARNING: Found equal dolic_c at edge jb, je=',jb, je, ' below dolic_e=', dol_e
    !              CALL message(TRIM(routine), TRIM(message_text))
    !            END IF
    !
    !          END IF ! iic1 > 0 .AND. iic2 > 0
    !
    !        END IF
    !
    !      END DO
    !    END DO
    !TODO review
    IF(MAXVAL(v_base%dolic_c)>n_zlev.OR.MINVAL(v_base%dolic_c)<0)THEN
      CALL message(TRIM(routine), TRIM('something wrong with dolic_c'))
      CALL finish(TRIM(routine),'something wrong with dolic_c')
    ENDIF
    IF(MAXVAL(v_base%dolic_e)>n_zlev.OR.MINVAL(v_base%dolic_e)<0)THEN
      CALL message(TRIM(routine), TRIM('something wrong with dolic_e'))
      CALL finish(TRIM(routine),'something wrong with dolic_e')
    ENDIF
    !-----------------------------
    ! real bathymetry should not be used since individual bottom layer thickness is not implemented
    ! set values of bathymetry to new non-individual dolic values

    ! DO jb = 1, alloc_cell_blocks
    !   i_endidx=nproma
    !   IF (jb==alloc_cell_blocks) i_endidx=npromz_c
    !     DO jc = 1, i_endidx
    !       v_base%bathymetry_c(jc,jb) = &
    !         &  -v_base%zlev_i(v_base%dolic_c(jc,jb)+1)
    !   ENDDO
    ! ENDDO

    ! DO jb = 1, nblks_e
    !   i_endidx=nproma
    !   IF (jb==nblks_e) i_endidx=npromz_e
    !     DO je = 1, i_endidx
    !       v_base%bathymetry_e(je,jb) = &
    !         &  -v_base%zlev_i(v_base%dolic_e(je,jb)+1)
    !   ENDDO
    ! ENDDO

    !$ACC ENTER DATA COPYIN(v_base, v_base%lsm_c)

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_ho_base

  !-------------------------------------------------------------------------
  !
  !
  !>
  !! Initializes 2-dimensional definitions of basins for the ocean state variable.
  !!
  !! The 2-dimensional land-sea-mask is used to define basins for calculation of
  !! meridional overturning circulation (MOC) and to define areas of certain interest.
  !!
  !!
  !!  no-mpi parallelized
!<Optimize:inUse>
  SUBROUTINE init_ho_basins( patch_2d, v_base )

    TYPE(t_patch), TARGET, INTENT(inout)       :: patch_2d
    TYPE(t_hydro_ocean_base), INTENT(inout)    :: v_base

    REAL(wp) :: z_sync_c(nproma,patch_2d%alloc_cell_blocks)
    REAL(wp) :: z_sync_e(nproma,patch_2d%nblks_e)
    INTEGER :: ibase   (nproma,patch_2d%alloc_cell_blocks)
    INTEGER :: iarea   (nproma,patch_2d%alloc_cell_blocks)

    INTEGER :: i_startidx_c, i_endidx_c
    INTEGER :: jb, jc, jk, alloc_cell_blocks, npromz_c, i
    INTEGER :: no_cor, g_cor, no_glb, jiter, iter
    INTEGER :: n_idx(3), n_blk(3)
    REAL(wp) :: z60n, z30n, z30s, z85s, z10n, z100w
    REAL(wp) :: z_lat_deg, z_lon_deg
    REAL(wp) :: z_lon_pta, z_lon_ati, z_lon_itp, z_lon_ind, z_lon_nam, z_lon_med

    TYPE(t_subset_range), POINTER :: all_cells, cells_in_domain
    LOGICAL :: is_area_5, is_area_8

    TYPE(t_ocean_regions) :: ocean_regions
    TYPE(t_ocean_basins)  :: ocean_basins

    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_initialization:init_ho_basins'

    !-----------------------------------------------------------------------------
    all_cells => patch_2d%cells%ALL
    cells_in_domain => patch_2d%cells%in_domain
    CALL message (TRIM(routine), 'start')

    z_sync_c(:,:) = 0.0_wp
    z_sync_e(:,:) = 0.0_wp

    ! values for the blocking
    alloc_cell_blocks = patch_2d%alloc_cell_blocks
    !nblks_e = patch_2D%nblks_e
    !nblks_v = patch_2D%nblks_v
    npromz_c = patch_2d%npromz_c
    !npromz_e = patch_2D%npromz_e



    !-----------------------------
    ! Define borders of region:
    !  two problematic regions remain that can be accessed via space filling curves
    !   - Caribbian Sea is partly divided in Pacific/Atlantic (border are land points)
    !     iterative loop is used as first guess, see below
    !   - Indonesian Region is both in Pacific/Indian Ocean
    !     there is no clear land border since the Indonesian Throughflow(s) exist
    !     here a slanted geographic line can be implemented

    z60n     =  61.0_wp

    z30n     =  30.0_wp
    z30s     = -30.0_wp
    z85s     = -85.0_wp
    z_lon_pta = -70.0_wp   !  Pac/Atl - Drake Passage
    z_lon_ati =  25.0_wp   !  Atl/Ind - South Africa
    z_lon_itp = 115.0_wp   !  Ind/Pac - Australia-Indonesia
    z_lon_ind = 100.0_wp   !  Ind/Pac - Indonesia (north of Equator)
    z_lon_nam = -90.0_wp   !  Pac/Atl - North America
    z_lon_med =  50.0_wp   !  Atl/Ind - Mediterranean

    ! for Caribbean: 5N and 105W are crucial
    z10n     =    5.0_wp
    z100w    = -105.0_wp

    !-----------------------------
    ! Fill ocean areas:
    iarea(:,:) = -1
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)

      DO jc = i_startidx_c, i_endidx_c

        IF (v_base%lsm_c(jc,1,jb) >= boundary) CYCLE

        ! get lat/lon of actual cell
        z_lat_deg = rad2deg * patch_2d%cells%center(jc,jb)%lat
        z_lon_deg = rad2deg * patch_2d%cells%center(jc,jb)%lon
        IF (z_lon_deg >  180.0_wp) z_lon_deg = z_lon_deg-360.0_wp
        IF (z_lon_deg < -180.0_wp) z_lon_deg = z_lon_deg+360.0_wp

        ! Arctic Ocean: default
        iarea(jc,jb) = ocean_regions%arctic_ocean

        ! GIN Sea (not yet)

        ! Labrador Sea (not yet)

        ! 4 = North Atlantic
        IF (                                                           &
          & (z_lat_deg >= z30n      .AND. z_lat_deg < z60n)     .AND.  &
          & (z_lon_deg >= z_lon_nam .AND. z_lon_deg < z_lon_med)       &
          & ) iarea(jc,jb) = ocean_regions%north_atlantic

        ! 9 = North Pacific
        IF (                                                           &
          & (z_lat_deg >= z30n      .AND. z_lat_deg < z60n)     .AND.  &
          & (z_lon_deg >= z_lon_itp .OR.  z_lon_deg < z_lon_nam)       &
          & ) iarea(jc,jb) = ocean_regions%north_pacific

        ! 5 = Tropical Atlantic (without Caribbean - yet)
        IF (                                                           &
          & (z_lat_deg >= z30s      .AND. z_lat_deg < z30n)     .AND.  &
          & (z_lon_deg >= z_lon_pta .AND. z_lon_deg < z_lon_med)       &
          & ) iarea(jc,jb) = ocean_regions%tropical_atlantic

        ! 6 = Southern Ocean
        IF (z_lat_deg < z30s) iarea(jc,jb) = ocean_regions%southern_ocean

        ! 7 = Indian (including Indonesian Pacific - yet)
        IF (                                                           &
          & (z_lat_deg >= z30s      .AND. z_lat_deg < z30n)     .AND.  &
          & (z_lon_deg >= z_lon_ati .AND. z_lon_deg < z_lon_itp)       &
          & ) iarea(jc,jb) = ocean_regions%indian_ocean

        ! 8 = Tropical Pacific
        IF (                                                           &
          & (z_lat_deg >= z30s      .AND. z_lat_deg < z30n)     .AND.  &
          & (z_lon_deg >= z_lon_itp .OR.  z_lon_deg < z_lon_pta)       &
          & ) iarea(jc,jb) = ocean_regions%tropical_pacific

        ! Crucial Region: Caribbean undefined (-33)
        IF (                                                           &
          & (z_lat_deg >= z10n      .AND. z_lat_deg < z30n)     .AND.  &
          & (z_lon_deg >= z100w     .AND. z_lon_deg < z_lon_pta)       &
          & ) iarea(jc,jb) = ocean_regions%caribbean

        ! Land points
        IF (v_base%lsm_c(jc,1,jb) >= boundary) iarea(jc,jb) = ocean_regions%land

      END DO
    END DO
    !chekc if iarea is the same
    z_sync_c(:,:) =  REAL(iarea(:,:),wp)
    CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
    iarea(:,:) = INT(z_sync_c(:,:))

    !-----------------------------
    ! Fill crucial areas:

    ! Caribbean: border is land point
    !  - border is at least one line of land points (not valid for Indonesian)
    !  - no cell has Pacific and Atlantic neighbor
    !  - iterative procedure (time consuming in higher resolution)

    ! disable p_test_run since iterations will be different
    CALL disable_sync_checks()

    g_cor=0
    !   Do jiter=1,1000 ! should be adjusted to higher resolution
    DO jiter=1,100
      !   Do jiter=1,10
      !   Do jiter=1,1

      no_cor = 0
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          is_area_5 = .FALSE.
          is_area_8 = .FALSE.
          IF (iarea(jc,jb) == ocean_regions%caribbean) THEN
            DO i = 1, 3  !  no_of_edges
              ! coordinates of neighbouring cells
              n_blk(i) = patch_2d%cells%neighbor_blk(jc,jb,i)
              n_idx(i) = patch_2d%cells%neighbor_idx(jc,jb,i)

              IF (n_idx(i) > 0) THEN
                IF (iarea(n_idx(i),n_blk(i)) == ocean_regions%tropical_atlantic) THEN       ! NEIGHBOR IS tROPaTL
                  is_area_5 = .TRUE.
                  EXIT
                END IF
              END IF

            END DO

            IF (is_area_5) THEN
              iarea(jc,jb)=ocean_regions%tropical_atlantic
              no_cor =  no_cor + 1
            ENDIF
          END IF
        END DO
      END DO

      no_glb = global_sum_array(no_cor)
      no_cor = no_glb
      g_cor=g_cor+no_cor

      iter=jiter
      IF (no_cor == 0) EXIT

      WRITE(message_text,'(a,i4,a,i8)') 'Corrected atlantic Caribbean region - iter=', &
        & jiter,' no of cor:',no_cor
      CALL message(TRIM(routine), TRIM(message_text))

      ! do sync

      z_sync_c(:,:) =  REAL(iarea(:,:),wp)
      CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
      iarea(:,:) = INT(z_sync_c(:,:))

    END DO

    WRITE(message_text,'(a,i4,a,i8)') 'Corrected atlantic Caribbean region - iterations=', &
      & iter,' no of cor:',g_cor
    CALL message(TRIM(routine), TRIM(message_text))

!    !chekc if iarea is the same
!    z_sync_c(:,:) =  REAL(iarea(:,:),wp)
!    CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
!    iarea(:,:) = INT(z_sync_c(:,:))
    !-----------------------------

    DO jiter=1,100
      !   Do jiter=1,10
      !   Do jiter=1,1

      no_cor = 0
      DO jb = cells_in_domain%start_block, cells_in_domain%end_block
        CALL get_index_range(cells_in_domain, jb, i_startidx_c, i_endidx_c)
        DO jc = i_startidx_c, i_endidx_c
          is_area_5 = .FALSE.
          is_area_8 = .FALSE.
          IF (iarea(jc,jb) == ocean_regions%caribbean) THEN
            DO i = 1, 3  !  no_of_edges
              ! coordinates of neighbouring cells
              n_blk(i) = patch_2d%cells%neighbor_blk(jc,jb,i)
              n_idx(i) = patch_2d%cells%neighbor_idx(jc,jb,i)

              IF (n_idx(i) > 0) THEN
                IF (iarea(n_idx(i),n_blk(i)) == ocean_regions%tropical_pacific) THEN  ! NEIGHBOR IS tROPpAC
                  is_area_8 = .TRUE.
                  EXIT
                END IF
              END IF

            END DO

            IF (is_area_8) THEN
              iarea(jc,jb)=ocean_regions%tropical_pacific
              no_cor =  no_cor + 1
            ENDIF
          END IF

        END DO
      END DO

      no_glb = global_sum_array(no_cor)
      no_cor = no_glb
      g_cor=g_cor+no_cor

      iter=jiter
      IF (no_cor == 0) EXIT

      WRITE(message_text,'(a,i4,a,i8)') 'Corrected pacific Caribbean region - iter=', &
        & jiter,' no of cor:',no_cor
      CALL message(TRIM(routine), TRIM(message_text))

      ! do sync

      z_sync_c(:,:) =  REAL(iarea(:,:),wp)
      CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
      iarea(:,:) = INT(z_sync_c(:,:))

    END DO

    WRITE(message_text,'(a,i4,a,i8)') 'Corrected total Caribbean region - iterations=', &
      & iter,' no of cor:',g_cor
    CALL message(TRIM(routine), TRIM(message_text))

    !chekc if iarea is the same
    CALL enable_sync_checks()
    z_sync_c(:,:) =  REAL(iarea(:,:),wp)
    CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
    iarea(:,:) = INT(z_sync_c(:,:))
    !-----------------------------

    !-----------------------------
    ! Fill ocean basins using ocean areas:

    WHERE ( iarea(:,:) <= ocean_regions%tropical_atlantic )
      ibase(:,:) = ocean_basins%atlantic
      ELSEWHERE ( iarea(:,:) == ocean_regions%southern_ocean )
      ibase(:,:) = ocean_regions%southern_ocean
      ELSEWHERE ( iarea(:,:) == ocean_regions%indian_ocean )
      ibase(:,:) = ocean_regions%indian_ocean
      ELSEWHERE ( iarea(:,:) >= ocean_regions%tropical_pacific)
      ibase(:,:) = ocean_basins%pacific
    END WHERE

    WHERE ( iarea(:,:) == ocean_regions%land )
      ibase(:,:) = ocean_regions%land
    END WHERE

    v_base%basin_c(:,:) = ibase(:,:)
    v_base%regio_c(:,:) = iarea(:,:)

    !-----------------------------
    ! set wet_c and wet_e to 1 at sea points including boundaries

    ! cells
    WHERE ( v_base%lsm_c(:,:,:) <= sea_boundary )
      v_base%wet_c(:,:,:) = 1.0_wp
    END WHERE

    ! edges
    WHERE ( v_base%lsm_e(:,:,:) <= sea_boundary )
      v_base%wet_e(:,:,:) = 1.0_wp
    END WHERE

    ! #slo# for test:
    !v_base%wet_c(:,:,:) = real(v_base%lsm_c(:,:,:),wp)
    !v_base%wet_e(:,:,:) = real(v_base%lsm_e(:,:,:),wp)

    ! intermediate levels: same as wet_c
    !WHERE ( v_base%lsm_c(:,:,:) <= SEA_BOUNDARY )
    !  v_base%wet_i(:,:,:) = 1.0_wp
    !END WHERE

    ! synchronize all elements of v_base - not necessary
    z_sync_c(:,:) =  REAL(v_base%basin_c(:,:),wp)
    CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
    v_base%basin_c(:,:) = INT(z_sync_c(:,:))
    z_sync_c(:,:) =  REAL(v_base%regio_c(:,:),wp)
    CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
    v_base%regio_c(:,:) = INT(z_sync_c(:,:))

    DO jk = 1, n_zlev

      z_sync_c(:,:) =  v_base%wet_c(:,jk,:)
      CALL sync_patch_array(sync_c, patch_2d, z_sync_c(:,:))
      v_base%wet_c(:,jk,:) = z_sync_c(:,:)

      z_sync_e(:,:) =  v_base%wet_e(:,jk,:)
      CALL sync_patch_array(sync_e, patch_2d, z_sync_e(:,:))
      v_base%wet_e(:,jk,:) = z_sync_e(:,:)

    END DO

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_ho_basins
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !! Modifies the already calculated Coriolis force, if beta-, f-plane or the nonrotating case
  !! is selected in the namelist. The tangent plane is associated to the center of the basin that is
  !! specified in the namelist. An alternative would be to associate it to the nearest edge/vertex,
  !! but this is not implemented yet, and i expect it to have a minor effect.
  !!
  !! The Coriolis parameter is specified for edges (needed in RBF-discretization) and at vertices
  !! (needed in mimetic discreization).
  !! The land-sea masks are not taken into account here. This would require to extend the
  !! 2D-coriolis-structure to a 3D one
  !!
  !!
!<Optimize:inUse>
  SUBROUTINE init_coriolis_oce( patch_2D )
    !
    IMPLICIT NONE
    !
    !
    TYPE(t_patch), TARGET, INTENT(inout) :: patch_2D
    !
    INTEGER :: jb, je, jv
    INTEGER :: StartEdgeIndex, EndEdgeIndex
    INTEGER :: StartVertexIndex, EndVertexIndex
    TYPE(t_geographical_coordinates) :: gc1,gc2
    TYPE(t_cartesian_coordinates) :: xx1, xx2
    REAL(wp) :: z_y, coriolis_lat
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = ('mo_ocean_initialization:init_coriolis_oce')
    TYPE(t_subset_range), POINTER :: all_verts, all_edges, all_cells
    !-----------------------------------------------------------------------
    all_verts => patch_2D%verts%ALL
    all_edges => patch_2D%edges%ALL
    all_cells => patch_2D%cells%ALL

    CALL message (TRIM(routine), 'start')

    SELECT CASE (coriolis_type)

    CASE(beta_plane_coriolis)

      CALL message (TRIM(routine), 'BETA_PLANE_CORIOLIS: set to linear approximation')

      coriolis_lat = basin_center_lat * deg2rad
      gc1%lat = basin_center_lat* deg2rad - 0.5_wp*basin_height_deg*deg2rad
      gc1%lon = 0.0_wp
      xx1=gc2cc(gc1)

      DO jb = all_verts%start_block, all_verts%end_block
        CALL get_index_range(all_verts, jb, StartVertexIndex, EndVertexIndex)
        DO jv = StartVertexIndex, EndVertexIndex
          !z_y = grid_sphere_radius*(patch_2D%verts%vertex(jv,jb)%lat - coriolis_lat)
          gc2%lat = patch_2D%verts%vertex(jv,jb)%lat!*deg2rad
          gc2%lon = 0.0_wp
          xx2=gc2cc(gc2)
          z_y = grid_sphere_radius * arc_length(xx2,xx1)
          patch_2D%verts%f_v(jv,jb) = 2.0_wp * grid_angular_velocity * &
            & ( SIN(coriolis_lat) + (COS(coriolis_lat)/grid_sphere_radius)*z_y)
          !  write(*,*)'beta', jv,jb,z_beta_plane_vort,2.0_wp*grid_angular_velocity*sin(coriolis_lat),&
          !  &2.0_wp*grid_angular_velocity*((cos(coriolis_lat)/grid_sphere_radius)*z_y)
        END DO
      END DO

      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, StartEdgeIndex, EndEdgeIndex)
        DO je = StartEdgeIndex, EndEdgeIndex
          ! depends on basin_center_lat only - not dependent on center_lon, basin_width or height
          gc2%lat = patch_2D%edges%center(je,jb)%lat!*deg2rad
          gc2%lon = 0.0_wp
          xx2=gc2cc(gc2)
          z_y = grid_sphere_radius*arc_length(xx2,xx1)

          !z_y = patch_2D%edges%center(je,jb)%lat - coriolis_lat
          patch_2D%edges%f_e(je,jb) = 2.0_wp * grid_angular_velocity * &
            & ( SIN(coriolis_lat) + &
            & (COS(coriolis_lat)/grid_sphere_radius)*z_y)
        END DO
      END DO
    CASE(f_plane_coriolis)

      CALL message (TRIM(routine), 'F_PLANE_CORIOLIS: set to constant value')

      coriolis_lat =  coriolis_fplane_latitude* deg2rad

      patch_2D%cells%f_c  = 2.0_wp*grid_angular_velocity*SIN(coriolis_lat)
      patch_2D%edges%f_e  = 2.0_wp*grid_angular_velocity*SIN(coriolis_lat)
      patch_2D%verts%f_v  = 2.0_wp*grid_angular_velocity*SIN(coriolis_lat)

    CASE(zero_coriolis)

      CALL message (TRIM(routine), 'ZERO_CORIOLIS: set to zero')
      patch_2D%cells%f_c = 0.0_wp
      patch_2D%verts%f_v = 0.0_wp
      patch_2D%edges%f_e = 0.0_wp

    CASE(full_coriolis)

      DO jb = all_verts%start_block, all_verts%end_block
        CALL get_index_range(all_verts, jb, StartVertexIndex, EndVertexIndex)
        DO jv = StartVertexIndex, EndVertexIndex
          patch_2D%verts%f_v(jv,jb) = 2.0_wp * grid_angular_velocity * SIN(patch_2D%verts%vertex(jv,jb)%lat)
        END DO
      END DO
      DO jb = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, jb, StartEdgeIndex, EndEdgeIndex)
        DO je = StartEdgeIndex, EndEdgeIndex
          patch_2D%edges%f_e(je,jb) = 2.0_wp * grid_angular_velocity * SIN(patch_2D%edges%center(je,jb)%lat)
        END DO
      END DO
       DO jb = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, jb, StartEdgeIndex, EndEdgeIndex)
        DO je = StartEdgeIndex, EndEdgeIndex
          patch_2D%cells%f_c(je,jb) = 2.0_wp * grid_angular_velocity * SIN(patch_2D%cells%center(je,jb)%lat)
        END DO
      END DO

    END SELECT

    CALL message (TRIM(routine), 'end')

  END SUBROUTINE init_coriolis_oce




  !-------------------------------------------------------------------------
  !>
  !! Allocation of basic 3-dimensional patch structure. This sbr assumes that
  !! the 2D horizontal patch components is already initialized.
  !
  !
  !!

!<Optimize:inUse>
  SUBROUTINE init_patch_3d(patch_3d, p_ext_data, v_base)

    TYPE(t_patch_3d ),TARGET, INTENT(inout)    :: patch_3d
    TYPE(t_external_data),    INTENT(inout)    :: p_ext_data
    TYPE(t_hydro_ocean_base), INTENT(inout)    :: v_base
    ! local variables
    INTEGER :: ist
    !     INTEGER :: alloc_cell_blocks, nblks_e, nblks_v, n_zlvp, n_zlvm!, ie
    INTEGER :: je,jc,jb,jk
    INTEGER :: i_startidx_c, i_endidx_c,StartEdgeIndex, EndEdgeIndex
    CHARACTER(LEN=max_char_length), PARAMETER :: &
      & routine = 'mo_ocean_initialization:init_patch_3D'

    TYPE(t_patch), POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells
    TYPE(t_subset_range), POINTER :: all_edges
    TYPE(t_subset_range), POINTER :: owned_cells
    TYPE(t_subset_range), POINTER :: edges_in_domain
    TYPE(t_subset_range), POINTER :: owned_verts

    INTEGER :: vertex_block, vertex_index, start_index, end_index, end_level
    INTEGER :: edge_block, edge_index, neighbor
    INTEGER :: land_edges, sea_edges, boundary_edges
    REAL(wp), ALLOCATABLE :: z_sync_v(:,:)
    REAL(wp), PARAMETER :: z_fac_limitthick = 0.8_wp  !  limits additional thickness of bottom cell
    REAL(wp) :: global_ocean_volume
    REAL(wp) :: z_prism_center_dist_e
    REAL(wp) :: temp_thick
    INTEGER, POINTER :: dolic_c(:,:), dolic_e(:,:)
    INTEGER :: max_edge_level
    INTEGER :: cell1_idx, cell1_blk, cell2_idx, cell2_blk
    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_initialization:init_patch_3D"

    !-----------------------------------------------------------------------------
    CALL message (TRIM(routine), 'start')
    patch_2d    => patch_3d%p_patch_2d(1)
    all_cells   => patch_2d%cells%ALL
    all_edges   => patch_2d%edges%ALL
    owned_cells => patch_2d%cells%owned
    edges_in_domain => patch_2d%edges%in_domain
    owned_verts => patch_2d%verts%owned
    !-------------------------------------------------------------------------

    !CALL message(TRIM(routine), 'start to construct basic hydro ocean state')

    !Copy indormation from v_base

    patch_3d%p_patch_1d(1)%zlev_i             = v_base%zlev_i
    patch_3d%p_patch_1d(1)%zlev_m             = v_base%zlev_m
    patch_3d%p_patch_1d(1)%del_zlev_i         = v_base%del_zlev_i
    patch_3d%p_patch_1d(1)%del_zlev_m         = v_base%del_zlev_m
    patch_3d%p_patch_1d(1)%inv_del_zlev_m(:)  = 0.0_wp
    DO jk = 1,n_zlev
      IF (v_base%del_zlev_m(jk) > 0.0_wp) &
        & patch_3d%p_patch_1d(1)%inv_del_zlev_m(jk) = 1.0_wp / v_base%del_zlev_m(jk)
    ENDDO

    patch_3d%p_patch_1d(1)%n_zlev            = v_base%n_zlev
    patch_3d%p_patch_1d(1)%n_zlvp            = v_base%n_zlvp
    patch_3d%p_patch_1d(1)%n_zlvm            = v_base%n_zlvm

    patch_3d%wet_e                           = v_base%wet_e
    patch_3d%wet_c                           = v_base%wet_c

    patch_3d%lsm_e                           = v_base%lsm_e
    patch_3d%lsm_c                           = v_base%lsm_c

    patch_3d%surface_cell_sea_land_mask(:,:) = p_ext_data%oce%lsm_ctr_c(:,:)
    patch_3d%surface_edge_sea_land_mask(:,:) = p_ext_data%oce%lsm_ctr_e(:,:)

    patch_3d%basin_c               = v_base%basin_c
    patch_3d%regio_c               = v_base%regio_c

    patch_3d%p_patch_1d(1)%dolic_c = v_base%dolic_c
    patch_3d%p_patch_1d(1)%dolic_e = v_base%dolic_e
    dolic_c => patch_3d%p_patch_1d(1)%dolic_c
    dolic_e => patch_3d%p_patch_1d(1)%dolic_e

    !dolic arrys should never contain 1, only 0 for land or values > 1
    IF(iswm_oce/=1.AND.n_zlev>1)THEN
      WHERE (dolic_c(:,:) < min_dolic) dolic_c(:,:) = 0
      WHERE (dolic_e(:,:) < min_dolic) dolic_e(:,:) = 0
    ENDIF
    patch_3d%p_patch_1d(1)%inv_prism_thick_c(:,:,:)       = 0.0_wp
    patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(:,:,:) = 0.0_wp
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        end_level = dolic_c(jc,jb)
        DO jk=1, end_level
          patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(jc,jk,jb) = v_base%del_zlev_m(jk)
          patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb)          = v_base%del_zlev_m(jk)
          patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,jb)    = v_base%del_zlev_i(jk)
          patch_3d%p_patch_1d(1)%prism_volume(jc,jk,jb)           = &
            & patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb) * &
            & patch_3d%p_patch_2d(1)%cells%area(jc,jb)
          patch_3d%p_patch_1d(1)%depth_CellMiddle(jc,jk,jb)      = patch_3d%p_patch_1d(1)%zlev_m(jk)
          patch_3d%p_patch_1d(1)%depth_CellInterface(jc,jk,jb)   = patch_3d%p_patch_1d(1)%zlev_i(jk)

          IF (patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb) > 0.0_wp) THEN
            patch_3d%p_patch_1d(1)%inv_prism_thick_c(jc,jk,jb) = &
              & 1.0_wp/patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb)
            patch_3d%p_patch_1d(1)%invConstantPrismThickness(jc,jk,jb) = &
              patch_3d%p_patch_1d(1)%inv_prism_thick_c(jc,jk,jb)
          ENDIF
          IF (v_base%del_zlev_i(jk) > 0.0_wp)  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)= 1.0_wp/v_base%del_zlev_i(jk)

          patch_3d%p_patch_1d(1)%constantPrismCenters_Zdistance(jc,jk,jb) = &
            & patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,jb)
          patch_3d%p_patch_1d(1)%constantPrismCenters_invZdistance(jc,jk,jb) = &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)

        END DO
        IF (end_level > 0) THEN
          patch_3d%p_patch_1d(1)%depth_CellInterface(jc,end_level+1,jb)   = patch_3d%p_patch_1d(1)%zlev_i(end_level+1)
          patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,end_level+1,jb)   = &
            & patch_3d%p_patch_1d(1)%prism_thick_c(jc,end_level,jb) * 0.5_wp
          IF (patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,end_level+1,jb) > 0.0_wp)  &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,end_level+1,jb)= &
              &   1.0_wp / patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,end_level+1,jb)
          patch_3d%p_patch_1d(1)%constantPrismCenters_Zdistance(jc,end_level+1,jb) = &
            & patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,end_level+1,jb)
          patch_3d%p_patch_1d(1)%constantPrismCenters_invZdistance(jc,end_level+1,jb) = &
            & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,end_level+1,jb)
        ENDIF

        ! set bottom/columns values
        jk = dolic_c(jc,jb)
        IF (jk >= min_dolic) THEN

          ! Bottom and column thickness for horizontally constant prism thickness
          patch_3d%bottom_thick_c(jc,jb) = patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb)
          patch_3d%column_thick_c(jc,jb) = v_base%zlev_i(jk+1)   !  lower bound of cell is at dolic_c+1

          ! Preliminary partial cells conform with l_max_bottom=false only
          IF (l_partial_cells) THEN

            ! Partial cell ends at real bathymetry below upper boundary zlev_i(dolic)
            ! at most one dry cell as neighbor is allowed, therefore bathymetry can be much deeper than corrected dolic
            ! maximum thickness limited to an additional part of the thickness of the underlying cell

            !! This assumes that bathymetry cannot be exactly the same as zlev_i
            !! So we put in a condition ensuring the cell thickness does not
            !! become too small
            IF (jk < n_zlev) THEN
              temp_thick =             &
                & MIN(-p_ext_data%oce%bathymetry_c(jc,jb)-v_base%zlev_i(jk), &
                & v_base%del_zlev_m(jk)+z_fac_limitthick*v_base%del_zlev_m(jk+1))
            ELSE
              ! maximum thickness limited to a similar factor of the thickness of the current cell
              temp_thick =             &
                & MIN(-p_ext_data%oce%bathymetry_c(jc,jb)-v_base%zlev_i(jk), &
                & (1.0_wp+z_fac_limitthick)*v_base%del_zlev_m(jk))
            ENDIF

            !! Only if the new thickness is above a threshold so we allow it to
            !! change
            IF (temp_thick .GT. 2.0_wp) THEN
              patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb) = temp_thick
            END IF


            !  this is necessary update for flat surface array but leads to abort in height equation
            patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(jc,jk,jb) =      &
              & patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb)
            patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,jb) = 0.5_wp* &
              & (v_base%del_zlev_m(jk-1) + patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb))
            patch_3d%p_patch_1d(1)%inv_prism_thick_c(jc,jk,jb)      =      &
              & 1.0_wp/patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb)
            patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc,jk,jb)=      &
              & 1.0_wp/patch_3d%p_patch_1d(1)%prism_center_dist_c(jc,jk,jb)

            patch_3d%p_patch_1d(1)%constantPrismCenters_Zdistance(jc, jk, jb) = &
              & patch_3d%p_patch_1d(1)%prism_center_dist_c(jc, jk, jb)
            patch_3d%p_patch_1d(1)%constantPrismCenters_invZdistance(jc, jk, jb)= &
              & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc, jk, jb)

            !  write(0,*)'XXXX flat_sfc_c',jk,jc,jb,&
            !    &patch_3D%p_patch_1D(1)%prism_thick_flat_sfc_c(jc,jk,jb),patch_3D%p_patch_1D(1)%del_zlev_m(jk)

            patch_3d%p_patch_1d(1)%constantPrismCenters_Zdistance(jc, jk, jb) = &
              & patch_3d%p_patch_1d(1)%prism_center_dist_c(jc, jk, jb)
            patch_3d%p_patch_1d(1)%constantPrismCenters_invZdistance(jc, jk, jb)= &
              & patch_3d%p_patch_1d(1)%inv_prism_center_dist_c(jc, jk, jb)

            ! bottom and column thickness for solver and output
            ! bottom cell thickness at jk=dolic
            patch_3d%bottom_thick_c(jc,jb) = patch_3d%p_patch_1d(1)%prism_thick_c(jc,jk,jb)
            ! column cell thickness: add upper column without elevation
            patch_3d%column_thick_c(jc,jb) = v_base%zlev_i(jk) + patch_3d%bottom_thick_c(jc,jb)

          ENDIF ! l_partial_cells
        ENDIF ! MIN_DOLIC

      END DO
    END DO

    patch_3d%p_patch_1d(1)%inv_prism_thick_e(:,:,:)       = 0.0_wp
    patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(:,:,:) = 0.0_wp
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, StartEdgeIndex, EndEdgeIndex)
      DO je = StartEdgeIndex, EndEdgeIndex


        cell1_idx = patch_2d%edges%cell_idx(je, jb, 1)
        cell1_blk = patch_2d%edges%cell_blk(je, jb, 1)
        cell2_idx = patch_2d%edges%cell_idx(je, jb, 2)
        cell2_blk = patch_2d%edges%cell_blk(je, jb, 2)


        DO jk=1, dolic_e(je,jb)
          patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,jk,jb) = v_base%del_zlev_m(jk)
          patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)          = v_base%del_zlev_m(jk)

          IF (v_base%del_zlev_m(jk) > 0.0_wp) &
            patch_3d%p_patch_1d(1)%inv_prism_thick_e(je,jk,jb)      = 1.0_wp/v_base%del_zlev_m(jk)
          IF (v_base%del_zlev_i(jk) > 0.0_wp) &
            patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(je,jk,jb)= 1.0_wp/v_base%del_zlev_i(jk)
        END DO
        DO jk = dolic_e(je,jb) + 1, n_zlev
          patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,jk,jb) = 0.0_wp
          patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)          = 0.0_wp
        ENDDO

        ! bottom/columns values
        jk = v_base%dolic_e(je,jb)
        IF (jk >= min_dolic) THEN

          ! Bottom and column thickness for horizontally constant prism thickness
          patch_3d%bottom_thick_e(je,jb) = patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
          patch_3d%column_thick_e(je,jb) = v_base%zlev_i(jk+1)   !  lower bound is below dolic_e

          ! Preliminary partial cells conform with l_max_bottom=false only
          IF (l_partial_cells) THEN

            !! FIXME: edge bathymetry seems to have bugs
            !! This is a simpler approximation
            patch_3D%p_patch_1D(1)%prism_thick_e(je, jk, jb) = &
              MIN(patch_3D%p_patch_1D(1)%prism_thick_c(cell1_idx, jk, cell1_blk), &
                & patch_3D%p_patch_1D(1)%prism_thick_c(cell2_idx, jk, cell2_blk))

!            ! Partial cell ends at real bathymetry below upper boundary
!            ! zlev_i(dolic) at most one dry cell as neighbor is allowed,
!            ! therefore bathymetry can be much deeper than corrected dolic
!            ! maximum thickness limited to an additional part of the thickness
!            ! of the underlying cell
!            IF (jk < n_zlev) THEN
!              patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb) =             &
!                & MIN(-p_ext_data%oce%bathymetry_e(je,jb)-v_base%zlev_i(jk), &
!                & v_base%del_zlev_m(jk)+z_fac_limitthick*v_base%del_zlev_m(jk+1))
!            ELSEIF(jk >= n_zlev) THEN
!!              ! maximum thickness limited to a similar factor of the thickness of the current cell
!!              patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb) =             &
!!                & MIN(-p_ext_data%oce%bathymetry_e(je,jb)-v_base%zlev_i(jk), &
!!                & (1.0_wp+z_fac_limitthick)*v_base%del_zlev_m(jk))
!              patch_3D%p_patch_1D(1)%prism_thick_e(je, jk, jb) = &
!              MIN(patch_3D%p_patch_1D(1)%prism_thick_c(cell1_idx, jk, cell1_blk), &
!                patch_3D%p_patch_1D(1)%prism_thick_c(cell2_idx, jk, cell2_blk))
!            ENDIF

!            !  this is necessary update for flat surface array but leads to abort in height equation

            patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_e(je,jk,jb) =      &
              & patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
            patch_3d%p_patch_1d(1)%inv_prism_thick_e(je,jk,jb)      = &
              & 1.0_wp/patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
            ! no matching prism_center_dist_e ?
            !     patch_3D%p_patch_1D(1)%prism_center_dist_e(je,jk,jb) = 0.5_wp* &
            !       & (v_base%del_zlev_m(jk-1) + patch_3D%p_patch_1D(1)%prism_thick_e(je,jk,jb))
            z_prism_center_dist_e = 0.5_wp*                           &
              & (v_base%del_zlev_m(jk-1) + patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb))
            patch_3d%p_patch_1d(1)%inv_prism_center_dist_e(je,jk,jb)= &
              & 1.0_wp/z_prism_center_dist_e

            ! bottom and column thickness for solver and output
            ! bottom edge thickness at jk=dolic
            patch_3d%bottom_thick_e(je,jb) = patch_3d%p_patch_1d(1)%prism_thick_e(je,jk,jb)
            ! column edge thickness: add upper column without elevation
            patch_3d%column_thick_e(je,jb) = v_base%zlev_i(jk) + patch_3d%bottom_thick_e(je,jb)

          ENDIF ! l_partial_cells
        ENDIF ! jk>=MIN_DOLIC

      END DO
    END DO

    ! set halo values to zero in specific arrays for calculating global sum with respect to lsm
    ! cells
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        patch_3d%wet_halo_zero_c(jc,:,jb) = 0.0_wp
      END DO
    END DO
    DO jb = owned_cells%start_block, owned_cells%end_block
      CALL get_index_range(owned_cells, jb, i_startidx_c, i_endidx_c)
      DO jc = i_startidx_c, i_endidx_c
        patch_3d%wet_halo_zero_c(jc,:,jb) = patch_3d%wet_c(jc,:,jb)
      END DO
    END DO
    ! edges
    DO jb = all_edges%start_block, all_edges%end_block
      CALL get_index_range(all_edges, jb, StartEdgeIndex, EndEdgeIndex)
      DO je = StartEdgeIndex, EndEdgeIndex
        patch_3d%wet_halo_zero_e(je,:,jb) = 0.0_wp
      END DO
    END DO
    DO jb = edges_in_domain%start_block, edges_in_domain%end_block
      CALL get_index_range(edges_in_domain, jb, StartEdgeIndex, EndEdgeIndex)
      DO je = StartEdgeIndex, EndEdgeIndex
        patch_3d%wet_halo_zero_e(je,:,jb) = patch_3d%wet_e(je,:,jb)
      END DO
    END DO

    ! calculate ocean area and ocean volume with respect to land-sea-mask:
    !  - global 3-dim sum of volume is in ocean_volume(n_zlev+1)
    global_ocean_volume = 0.0_wp
    DO jk = 1, n_zlev
      patch_3d%p_patch_1d(1)%ocean_area(jk)   = global_sum_array( &
        & patch_2d%cells%area(:,:) * &
        & patch_3d%wet_halo_zero_c(:,jk,:) )
      patch_3d%p_patch_1d(1)%ocean_volume(jk) = global_sum_array( &
        & patch_2d%cells%area(:,:) * &
        & patch_3d%wet_halo_zero_c(:,jk,:)         * &
        & patch_3d%p_patch_1d(1)%prism_thick_flat_sfc_c(:,jk,:) )
      global_ocean_volume = global_ocean_volume + patch_3d%p_patch_1d(1)%ocean_volume(jk)
    END DO
    !-------------------------------------------------
    patch_3d%p_patch_1d(1)%ocean_volume(n_zlev+1) = global_ocean_volume

    !-------------------------------------------------
   ! calculate surface_vertex_sea_land_mask and vertex_bottomLevel
    DO vertex_block = owned_verts%start_block, owned_verts%end_block
      CALL get_index_range(owned_verts, vertex_block, start_index, end_index)
      DO vertex_index = start_index, end_index
        land_edges     = 0
        sea_edges      = 0
        boundary_edges = 0
        max_edge_level = 0

        DO neighbor=1, patch_2d%verts%num_edges(vertex_index,vertex_block)
          edge_index = patch_2d%verts%edge_idx(vertex_index, vertex_block, neighbor)
          edge_block = patch_2d%verts%edge_blk(vertex_index, vertex_block, neighbor)

          IF (edge_index > 0) THEN ! this should not be necessary

            SELECT CASE(patch_3d%lsm_e(edge_index, 1, edge_block))
            CASE(sea, sea_boundary)
              sea_edges = sea_edges + 1
            CASE(boundary)
              boundary_edges = boundary_edges + 1
            CASE(land, land_boundary)
              land_edges = land_edges + 1
            CASE default
              CALL finish(routine, "Uknown patch_3D%lsm_e" )
            END SELECT

            max_edge_level = MAX(max_edge_level, dolic_e(edge_index, edge_block))

          ENDIF

        ENDDO ! neighbor

        !        This is not true when land points are missing
        !        IF( MOD(boundary_edges,2) /= 0 ) THEN
        !          CALL finish (method_name,'MOD(boundary_edges,2) /= 0 !!')
        !        ENDIF

        patch_3d%surface_vertex_sea_land_mask(vertex_index, vertex_block)   = sea
        IF (boundary_edges > 0) THEN
          patch_3d%surface_vertex_sea_land_mask(vertex_index, vertex_block) = boundary
        ELSEIF (land_edges > 0) THEN
          patch_3d%surface_vertex_sea_land_mask(vertex_index, vertex_block) = land
          ! consistency check
          IF (sea_edges > 0) &
            & CALL finish(routine, "Inconsistent patch_3D%lsm_e" )
        ENDIF

        ! the vertex number of levels is the max levels of intersecting edges
        patch_3d%p_patch_1d(1)%vertex_bottomLevel(vertex_index, vertex_block) = max_edge_level

      ENDDO ! vertex_index
    ENDDO ! vertex_block
    ! sync the results
    ALLOCATE(z_sync_v(nproma,patch_2d%nblks_v),stat=ist)
    IF (ist /= success) THEN
      CALL finish (routine,'allocating surface_vertex_sea_land_mask failed')
    ENDIF
    z_sync_v(:,:) =  REAL(patch_3d%surface_vertex_sea_land_mask(:,:),wp)
    CALL sync_patch_array(sync_v, patch_2d, z_sync_v(:,:))
    patch_3d%surface_vertex_sea_land_mask(:,:) = INT(z_sync_v(:,:))
    z_sync_v(:,:) =  REAL(patch_3d%p_patch_1d(1)%vertex_bottomLevel(:,:),wp)
    CALL sync_patch_array(sync_v, patch_2d, z_sync_v(:,:))
    patch_3d%p_patch_1d(1)%vertex_bottomLevel(:,:) = INT(z_sync_v(:,:))
    DEALLOCATE(z_sync_v)
    !---------------------------------------

    CALL complete_ocean_subsets(patch_3d)

    CALL ocean_check_level_sea_land_mask(patch_3d)

   CALL dbg_print('init_patch_3d:thick_e',patch_3D%p_patch_1d(1)%prism_thick_e,&
      & "",4, in_subset=edges_in_domain)


  END SUBROUTINE init_patch_3d
  !------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE complete_ocean_subsets(patch_3d)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d

    CALL set_subset_ocean_vertical_layers(patch_3d)
    CALL ocean_subsets_ignore_land(patch_3d)
    CALL check_ocean_subsets(patch_3d)

  END SUBROUTINE complete_ocean_subsets
  !------------------------------------------------------------------------------------


  !------------------------------------------------------------------------------------
!<Optimize:inUse>
  SUBROUTINE set_subset_ocean_vertical_layers(patch_3d)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d

    TYPE(t_grid_cells),   POINTER :: cells
    TYPE(t_grid_edges),   POINTER :: edges
    INTEGER, POINTER :: cells_vertical_levels(:,:), edges_vertical_levels(:,:)

    INTEGER :: BLOCK, startidx, endidx, idx
    !-----------------------------------------------------------------------------

    cells  => patch_3d%p_patch_2d(1)%cells
    edges  => patch_3d%p_patch_2d(1)%edges

    cells_vertical_levels => patch_3d%p_patch_1d(1)%dolic_c
    edges_vertical_levels => patch_3d%p_patch_1d(1)%dolic_e

    cells%ALL%vertical_levels                => cells_vertical_levels
    cells%owned%vertical_levels              => cells_vertical_levels
    cells%in_domain%vertical_levels          => cells_vertical_levels
    cells%not_owned%vertical_levels          => cells_vertical_levels
    cells%not_in_domain%vertical_levels      => cells_vertical_levels
    cells%one_edge_in_domain%vertical_levels => cells_vertical_levels

    edges%ALL%vertical_levels           => edges_vertical_levels
    edges%in_domain%vertical_levels     => edges_vertical_levels
    edges%owned%vertical_levels         => edges_vertical_levels
    edges%not_owned%vertical_levels     => edges_vertical_levels
    edges%not_in_domain%vertical_levels => edges_vertical_levels

    cells%ALL%max_vertical_levels                = n_zlev
    cells%owned%max_vertical_levels              = n_zlev
    cells%in_domain%max_vertical_levels          = n_zlev
    cells%not_owned%max_vertical_levels          = n_zlev
    cells%not_in_domain%max_vertical_levels      = n_zlev
    cells%one_edge_in_domain%max_vertical_levels = n_zlev

    edges%ALL%max_vertical_levels           = n_zlev
    edges%in_domain%max_vertical_levels     = n_zlev
    edges%owned%max_vertical_levels         = n_zlev
    edges%not_owned%max_vertical_levels     = n_zlev
    edges%not_in_domain%max_vertical_levels = n_zlev

  END SUBROUTINE set_subset_ocean_vertical_layers
  !------------------------------------------------------------------------------------

  !------------------------------------------------------------------------------------
  ! ignore land points in range subsets
  ! this only works if the land points are re-ordered
!<Optimize:inUse>
  SUBROUTINE ocean_subsets_ignore_land(patch_3d)
    TYPE(t_patch_3d ),TARGET, INTENT(inout) :: patch_3d

    TYPE(t_patch),     POINTER :: patch_2d
    TYPE(t_subset_range), POINTER :: all_cells, all_edges, all_verts

    INTEGER :: BLOCK, startidx, endidx, idx
    !-----------------------------------------------------------------------------

    patch_2d => patch_3d%p_patch_2d(1)
    all_cells   => patch_2d%cells%ALL
    all_edges   => patch_2d%edges%ALL
    all_verts   => patch_2d%verts%ALL

    !--------------------------------------------------------------------------------
    ! exclude land from subsets, if requested
    IF ( .FALSE.) THEN

      ! cells
      DO BLOCK = all_cells%start_block, all_cells%end_block
        CALL get_index_range(all_cells, BLOCK, startidx, endidx)
        DO idx = startidx, endidx
          IF (patch_3d%surface_cell_sea_land_mask(idx, BLOCK) >= 0) THEN
            patch_2d%cells%decomp_info%halo_level(idx, BLOCK) = -1 !halo_levels_ceiling+1
            WRITE(0,*) "Removed land cell at ", idx, BLOCK
          ENDIF
        ENDDO
      ENDDO
      ! recalculate cells subsets
      CALL fill_subset(subset=patch_2D%cells%all, patch=patch_2D, &
        & mask=patch_2D%cells%decomp_info%halo_level, start_mask=0, end_mask=halo_levels_ceiling, located=on_cells)
      CALL fill_subset(subset=patch_2D%cells%owned, patch=patch_2D, &
        & mask=patch_2D%cells%decomp_info%halo_level, start_mask=0, end_mask=0, located=on_cells)

      CALL fill_subset(subset=patch_2D%cells%in_domain, patch=patch_2D, &
        & mask=patch_2D%cells%decomp_info%halo_level, start_mask=0, end_mask=0, located=on_cells)

      ! edges
      DO BLOCK = all_edges%start_block, all_edges%end_block
        CALL get_index_range(all_edges, BLOCK, startidx, endidx)
        DO idx = startidx, endidx
          IF (patch_3d%surface_edge_sea_land_mask(idx, BLOCK) == land) THEN
            patch_2d%edges%decomp_info%halo_level(idx, BLOCK) = -1 !halo_levels_ceiling+1
            WRITE(0,*) "Removed land edge at ", idx, BLOCK
          ENDIF
        ENDDO
      ENDDO
      ! recalculate edges subsets
      CALL fill_subset(subset=patch_2D%edges%all, patch=patch_2D, &
        & mask=patch_2D%edges%decomp_info%halo_level, start_mask=0, end_mask=halo_levels_ceiling, located=on_edges)
      CALL fill_subset(subset=patch_2D%edges%owned, patch=patch_2D, &
        & mask=patch_2D%edges%decomp_info%halo_level, start_mask=0, end_mask=0, located=on_edges)
      CALL fill_subset(subset=patch_2D%edges%in_domain, patch=patch_2D, &
        & mask=patch_2D%edges%decomp_info%halo_level, start_mask=0, end_mask=1, located=on_edges)

      ! verts
      DO BLOCK = all_verts%start_block, all_verts%end_block
        CALL get_index_range(all_verts, BLOCK, startidx, endidx)
        DO idx = startidx, endidx
          IF (patch_3d%surface_vertex_sea_land_mask(idx, BLOCK) == land) THEN
            patch_2d%verts%decomp_info%halo_level(idx, BLOCK) = -1 !halo_levels_ceiling+1
            WRITE(0,*) "Removed land vert at ", idx, BLOCK
          ENDIF
        ENDDO
      ENDDO
      ! recalculate verts subsets
      CALL fill_subset(subset=patch_2D%verts%all,  patch=patch_2D, &
        & mask=patch_2D%verts%decomp_info%halo_level, start_mask=0, end_mask=halo_levels_ceiling, located=on_vertices)
      CALL fill_subset(subset=patch_2D%verts%owned, patch=patch_2D, &
        & mask=patch_2D%verts%decomp_info%halo_level, start_mask=0, end_mask=0, located=on_vertices)
      CALL fill_subset(subset=patch_2D%verts%in_domain, patch=patch_2D, &
        & mask=patch_2D%verts%decomp_info%halo_level, start_mask=0, end_mask=1, located=on_vertices)

    ENDIF
  END SUBROUTINE ocean_subsets_ignore_land
  !------------------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
!<Optimize:inUse>
  SUBROUTINE set_del_zlev(n_zlev, dzlev_m, del_zlev_i, del_zlev_m, zlev_i, zlev_m)
    INTEGER,  INTENT(IN) :: n_zlev
    REAL(wp), INTENT(IN) :: dzlev_m(max_allocated_levels)
    REAL(wp)             :: del_zlev_i(n_zlev), del_zlev_m(n_zlev)
    REAL(wp)             :: zlev_i(n_zlev+1)  , zlev_m(n_zlev)

    INTEGER :: jk
  !!-------------------------------------
    CALL set_zlev(zlev_i, zlev_m, n_zlev, dzlev_m)
    ! del_zlev_i: distance between two z-coordinate surfaces.
    !             The first is the distance from the ocean surface = zlev_m(1)
    del_zlev_i(1) = zlev_m(1)
    DO jk = 2, n_zlev
      del_zlev_i(jk) = zlev_m(jk) -  zlev_m(jk-1)
    END DO

    del_zlev_m(:) = dzlev_m(1:n_zlev)

  END SUBROUTINE set_del_zlev
  !------------------------------------------------------------------------------------


!<Optimize:inUse>
  FUNCTION is_initial_timestep(timestep)
    INTEGER :: timestep
    LOGICAL is_initial_timestep

    IF (timestep == 1 .AND. .NOT. isRestart()) THEN
      is_initial_timestep = .TRUE.
    ELSE
      is_initial_timestep = .FALSE.
    END IF
  END FUNCTION is_initial_timestep

  !----------------------------------------------------------------------------

END MODULE mo_ocean_initialization
