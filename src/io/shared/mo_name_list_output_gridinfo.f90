!
! Module handling the collection of grid information for synchronous
! and asynchronous output.
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

MODULE mo_name_list_output_gridinfo

  USE mo_cdi,                               ONLY: DATATYPE_PACK16, TIME_CONSTANT, TSTEP_CONSTANT, vlistDefVar, cdiEncodeParam, &
                                                & streamWriteVar, vlistDefVarDatatype, vlistDefVarName, vlistDefVarTsteptype, &
                                                & vlistDefVarParam, gridDefXvals, gridDefYvals, gridDefXbounds, gridDefYbounds, &
                                                & GRID_UNSTRUCTURED
  USE mo_zaxis_type,                        ONLY: ZA_surface
  USE mo_kind,                              ONLY: wp
  USE mo_parallel_config,                   ONLY: nproma
  USE mo_exception,                         ONLY: finish
  USE mo_model_domain,                      ONLY: t_patch
  USE mo_math_types,                        ONLY: t_geographical_coordinates
  USE mo_math_utilities,                    ONLY: check_orientation
  USE mo_communication,                     ONLY: t_comm_gather_pattern, exchange_data,     &
    &                                             t_comm_allgather_pattern,       &
    &                                             setup_comm_allgather_pattern,   &
    &                                             delete_comm_allgather_pattern,   &
    &                                             setup_comm_gather_pattern,   &
    &                                             delete_comm_gather_pattern
  USE mo_grib2,                             ONLY: t_grib2_var, grib2_var
  USE mo_grib2_util,                        ONLY: set_GRIB2_additional_keys,                &
    &                                             set_GRIB2_ensemble_keys,                  &
    &                                             set_GRIB2_local_keys
  USE mo_lonlat_grid,                       ONLY: t_lon_lat_grid, compute_lonlat_specs,     &
    &                                             rotate_latlon_grid
  USE mo_intp_lonlat_types,                 ONLY: lonlat_grids
  USE mo_math_constants,                    ONLY: pi, rad2deg
  USE mo_impl_constants,                    ONLY: SUCCESS, min_rlcell_int, min_rledge_int,  &
    &                                             min_rlvert, vname_len
  USE mo_cdi_constants,                     ONLY: GRID_CELL
  USE mo_mpi,                               ONLY: p_comm_work_2_io,                         &
    &                                             my_process_is_io, &
    &                                             my_process_is_mpi_workroot
  USE mo_master_control,                    ONLY: my_process_is_oceanic
  USE mo_gribout_config,                    ONLY: gribout_config
  USE mo_loopindices,                       ONLY: get_indices_c, get_indices_e, get_indices_v
  USE mo_util_string,                       ONLY: one_of
  USE mo_name_list_output_types,            ONLY: t_patch_info, t_grid_info, t_output_file, &
    &                                             REMAP_NONE, REMAP_REGULAR_LATLON,         &
    &                                             ILATLON, ICELL, IEDGE, IVERT, IRLAT,      &
    &                                             IRLON, GRB2_GRID_INFO_NAME
  USE mo_name_list_output_zaxes_types,      ONLY: t_verticalAxis
  USE mo_var, ONLY: level_type_ml, level_type_pl, level_type_hl, level_type_il
  USE mo_reorder_info,                      ONLY: t_reorder_info
  USE mo_netcdf
#ifdef HAVE_CDI_PIO
  USE mo_reorder_info,                      ONLY: ri_cpy_blk2part
  USE mo_parallel_config,                   ONLY: pio_type
  USE mo_impl_constants,                    ONLY: pio_type_cdipio
  USE yaxt,                                 ONLY: xt_idxlist, xt_stripe, &
    &                                             xt_idxstripes_new
  USE ppm_extents,                          ONLY: extent
  USE mo_decomposition_tools,               ONLY: uniform_partition_start, uniform_partition
  USE mo_mpi,                               ONLY: p_n_work, p_pe_work
  USE iso_c_binding,                        ONLY: c_int
#endif
  IMPLICIT NONE

  PRIVATE

  !-----------------------------------------------------------------
  ! include NetCDF headers (direct NetCDF library calls are required
  ! for output of grid information).
  !-----------------------------------------------------------------

  ! constants
  PUBLIC :: GRID_INFO_NONE, GRID_INFO_FILE, GRID_INFO_BCAST
  ! subroutines
  PUBLIC :: deallocate_all_grid_info
  PUBLIC :: collect_all_grid_info
#ifdef HAVE_CDI_PIO
  PUBLIC :: distribute_all_grid_info
#endif
  PUBLIC :: set_grid_info_netcdf
  PUBLIC :: set_grid_info_grb2
  PUBLIC :: copy_grid_info
  PUBLIC :: allgather_grid_info
  PUBLIC :: write_grid_info_grb2

  ! constants: modes how to collect grid information (for output)
  INTEGER, PARAMETER :: GRID_INFO_NONE  = 1
  INTEGER, PARAMETER :: GRID_INFO_FILE  = 2
  INTEGER, PARAMETER :: GRID_INFO_BCAST = 3


  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_name_list_output_gridinfo'

CONTAINS

  SUBROUTINE deallocate_grid_info(grid_info)
    TYPE(t_grid_info),  INTENT(INOUT) :: grid_info
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_grid_info"
    ! local variables
    INTEGER :: ierrstat

    IF (ALLOCATED(grid_info%lon)) THEN
      DEALLOCATE(grid_info%lon,  grid_info%lat,  &
        &        grid_info%lonv, grid_info%latv, &
        &        STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
    IF (ALLOCATED(grid_info%log_dom_starts)) THEN
      DEALLOCATE(grid_info%log_dom_starts, grid_info%log_dom_counts, &
        &        STAT=ierrstat)
      IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
    END IF
  END SUBROUTINE deallocate_grid_info


  SUBROUTINE deallocate_all_grid_info(patch_info)
    TYPE(t_patch_info),   INTENT(INOUT) :: patch_info
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::deallocate_grid_info"

    INTEGER :: i

    DO i = 1, 3
      CALL deallocate_grid_info(patch_info%grid_info(i))
    END DO
  END SUBROUTINE deallocate_all_grid_info


  SUBROUTINE collect_all_grid_info(p_patch, patch_info)
    TYPE(t_patch),        INTENT(IN)    :: p_patch
    TYPE(t_patch_info),   INTENT(INOUT) :: patch_info
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::collect_all_grid_info"
    INTEGER                             :: ierrstat, max_cell_connectivity, max_vertex_connectivity
    REAL(wp), ALLOCATABLE               :: lonv(:,:,:), latv(:,:,:)

    ! logical domain ID
    max_cell_connectivity   = p_patch%cells%max_connectivity
    max_vertex_connectivity = p_patch%verts%max_connectivity
    !-- collect domain data on working PE 0
    ! --cells
    ALLOCATE(lonv(nproma, p_patch%nblks_c, max_cell_connectivity), &
      &      latv(nproma, p_patch%nblks_c, max_cell_connectivity), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
!$omp parallel
    CALL cf_1_1_grid_cells(p_patch, lonv, latv)
!$omp end parallel
    CALL collect_grid_info(patch_info%nblks_glb_c,              &
      &                    p_patch%nblks_c,                     &
      &                    p_patch%cells%center,                &
      &                    lonv, latv,                          &
      &                    patch_info%grid_info(icell),         &
      &                    max_cell_connectivity,               &
      &                    patch_info%p_pat_c)
    DEALLOCATE(lonv, latv, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    !-- edges
    ALLOCATE(lonv(nproma, p_patch%nblks_e, 4), &
      &      latv(nproma, p_patch%nblks_e, 4), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
!$omp parallel
    CALL cf_1_1_grid_edges(p_patch, lonv, latv)
!$omp end parallel
    CALL collect_grid_info(patch_info%nblks_glb_e,              &
      &                    p_patch%nblks_e,                     &
      &                    p_patch%edges%center,                &
      &                    lonv, latv,                          &
      &                    patch_info%grid_info(iedge),         &
      &                    4,                                   &
      &                    patch_info%p_pat_e)
    DEALLOCATE(lonv, latv, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

    !-- verts
    ALLOCATE(lonv(nproma, p_patch%nblks_v, max_vertex_connectivity), &
      &      latv(nproma, p_patch%nblks_v, max_vertex_connectivity), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

!$omp parallel
    IF (my_process_is_oceanic()) THEN
      CALL cf_1_1_grid_verts_ocean(p_patch, lonv, latv)
    ELSE
      CALL cf_1_1_grid_verts(p_patch, lonv, latv)
    ENDIF
!$omp end parallel

    CALL collect_grid_info(patch_info%nblks_glb_v,              &
      &                    p_patch%nblks_v,                     &
      &                    p_patch%verts%vertex,                &
      &                    lonv, latv,                          &
      &                    patch_info%grid_info(ivert),          &
      &                    max_vertex_connectivity,             &
      &                    patch_info%p_pat_v)
    DEALLOCATE(lonv, latv, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')
  END SUBROUTINE collect_all_grid_info

#ifdef HAVE_CDI_PIO
  SUBROUTINE distribute_all_grid_info(p_patch, patch_info)
    TYPE(t_patch), TARGET, INTENT(in)   :: p_patch
    TYPE(t_patch_info),   INTENT(inout) :: patch_info

    CHARACTER(LEN=*), PARAMETER :: &
      routine = modname//"::distributed_all_grid_info"
    INTEGER :: ierrstat, max_a_size, igeom
    REAL(wp), ALLOCATABLE :: lonv(:), latv(:)
    TYPE cf_1_1_grid_ptr
      PROCEDURE(cf_1_1_grid_verts), NOPASS, POINTER :: p
    END TYPE cf_1_1_grid_ptr
    TYPE(cf_1_1_grid_ptr) :: cf_1_1(3)
    TYPE p_geographical_coordinates
      TYPE(t_geographical_coordinates), POINTER :: p(:,:)
    END TYPE p_geographical_coordinates
    TYPE(p_geographical_coordinates) :: p_coords(3)
    INTEGER :: max_conn(3), nblks(3)

    cf_1_1(icell)%p => cf_1_1_grid_cells
    cf_1_1(iedge)%p => cf_1_1_grid_edges
    IF (my_process_is_oceanic()) THEN
      cf_1_1(ivert)%p => cf_1_1_grid_verts_ocean
    ELSE
      cf_1_1(ivert)%p => cf_1_1_grid_verts
    ENDIF

    p_coords(icell)%p => p_patch%cells%center
    p_coords(iedge)%p => p_patch%edges%center
    p_coords(ivert)%p => p_patch%verts%vertex

    max_conn(icell) = p_patch%cells%max_connectivity
    max_conn(iedge) = 4
    max_conn(ivert) = p_patch%verts%max_connectivity

    nblks(icell) = p_patch%nblks_c
    nblks(iedge) = p_patch%nblks_e
    nblks(ivert) = p_patch%nblks_v

    DO igeom = icell, ivert
      CALL alloc_distributed_grid_info(patch_info%grid_info(igeom), &
        patch_info%ri(igeom)%n_own, max_conn(igeom))
    END DO

    max_a_size = nproma * MAXVAL(nblks * max_conn)
    ALLOCATE(lonv(max_a_size), latv(max_a_size))
!$omp parallel private(igeom)
    DO igeom = icell, ivert
      CALL create_distributed_grid_info(p_patch, patch_info%ri(igeom), &
        &                               lonv, latv, p_coords(igeom)%p, &
        &                               patch_info%grid_info(igeom), &
        &                               max_conn(igeom), nblks(igeom), &
        &                               cf_1_1(igeom)%p)
    END DO
!$omp end parallel
  END SUBROUTINE distribute_all_grid_info

  SUBROUTINE alloc_distributed_grid_info(grid_info, n_own, nconn)
    TYPE(t_grid_info), INTENT(inout) :: grid_info
    INTEGER, INTENT(in) :: n_own, nconn

    CHARACTER(len=*), PARAMETER :: &
         routine = modname//'::alloc_distributed_grid_info'
    INTEGER :: ierrstat

    ALLOCATE(grid_info%lon(n_own), grid_info%lat(n_own), &
      &      grid_info%lonv(nconn, n_own), grid_info%latv(nconn, n_own), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')
  END SUBROUTINE alloc_distributed_grid_info

  SUBROUTINE create_distributed_grid_info(p_patch, ri, lonv, latv, coordinates,&
    &                                     grid_info, nconn, nblks, cf_1_1_grid)
    TYPE(t_patch), INTENT(in) :: p_patch
    TYPE(t_reorder_info), INTENT(in) :: ri
    REAL(wp), TARGET, INTENT(inout) :: lonv(:), latv(:)
    TYPE(t_geographical_coordinates), INTENT(IN) :: coordinates(:,:)
    TYPE(t_grid_info), INTENT(inout) :: grid_info
    INTEGER, INTENT(in) :: nconn, nblks
    INTERFACE
      SUBROUTINE cf_1_1_grid(p_patch, lonv, latv)
        USE mo_model_domain, ONLY: t_patch
        USE mo_kind,         ONLY: wp
        TYPE(t_patch),      INTENT(IN)    :: p_patch
        REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
      END SUBROUTINE cf_1_1_grid
    END INTERFACE

    CHARACTER(len=*), PARAMETER :: &
      routine = modname//'::create_distributed_grid_info'
    INTEGER :: ierrstat, j, n_own, ofs
    REAL(wp), POINTER :: plonv(:,:,:), platv(:,:,:)

    n_own = ri%n_own
    plonv(1:nproma, 1:nblks, 1:nconn) => lonv
    platv(1:nproma, 1:nblks, 1:nconn) => latv
    CALL cf_1_1_grid(p_patch, plonv, platv)
!$omp barrier
    ofs = 0 ; CALL ri_cpy_blk2part(ri, coordinates%lon, grid_info%lon, ofs)
    ofs = 0 ; CALL ri_cpy_blk2part(ri, coordinates%lat, grid_info%lat, ofs)
    DO j = 1, nconn
      ofs = 0 ; CALL ri_cpy_blk2part(ri, plonv(:,:,j), grid_info%lonv(j,:), ofs)
      ofs = 0 ; CALL ri_cpy_blk2part(ri, platv(:,:,j), grid_info%latv(j,:), ofs)
    END DO
!$omp barrier
  END SUBROUTINE create_distributed_grid_info
#endif

  !------------------------------------------------------------------------------------------------
  !> SUBROUTINE collect_grid_info
  !
  !  Prepare the output of grid information: For each physical/logical
  !  patch we must collect the geographical locations of cells, edges,
  !  and vertices first on working PE 0 - from where it will be
  !  broadcast to the pure I/O PEs.
  !
  SUBROUTINE collect_grid_info(nblks_glb, nblks_loc, in_lonlat, lonv, latv, out_lonlat, &
    &                          dim3, p_pat)
    INTEGER,                          INTENT(IN)    :: nblks_glb,     &  ! global number of blocks
      &                                                nblks_loc         ! local  number of blocks
    TYPE(t_geographical_coordinates), INTENT(IN)    :: in_lonlat(:,:)
    REAL(wp),                         INTENT(IN)    :: lonv(:,:,:), latv(:,:,:)
    TYPE(t_grid_info),                INTENT(INOUT) :: out_lonlat
    INTEGER,                          INTENT(IN)    :: dim3
    TYPE(t_comm_gather_pattern),  POINTER           :: p_pat
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine  =  modname//"::collect_grid_info"
    INTEGER               :: ierrstat, jb, jc, idim
    REAL(wp), ALLOCATABLE :: r_tmp_lon(:,:), r_tmp_lat(:,:)

    ! allocate destination (on work root)
    IF ( my_process_is_mpi_workroot() ) THEN
      ALLOCATE(out_lonlat%lon (nproma*nblks_glb),      &
        &      out_lonlat%lat (nproma*nblks_glb),      &
        &      out_lonlat%lonv(dim3,nproma*nblks_glb), &
        &      out_lonlat%latv(dim3,nproma*nblks_glb), &
        &      STAT=ierrstat)
      out_lonlat%lon = 0._wp
      out_lonlat%lat = 0._wp
      out_lonlat%lonv = 0._wp
      out_lonlat%latv = 0._wp
    ELSE
      ALLOCATE(out_lonlat%lon(1), out_lonlat%lat (1),            &
        &      out_lonlat%lonv(dim3,1), out_lonlat%latv(dim3,1), &
        &      STAT=ierrstat)
    END IF
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    ! allocate temporary fields:
    ALLOCATE(r_tmp_lon (nproma, nblks_loc), r_tmp_lat (nproma, nblks_loc), &
      &      STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'ALLOCATE failed.')

    !-- part 1: exchange lon/lat coordinates:
    ! copy coordinates into sender array:
    DO jb=1,nblks_loc
      DO jc=1,nproma
        r_tmp_lon(jc,jb) = in_lonlat(jc,jb)%lon
        r_tmp_lat(jc,jb) = in_lonlat(jc,jb)%lat
      END DO
    END DO

    CALL exchange_data(in_array=r_tmp_lon(:,:), out_array=out_lonlat%lon(:), &
      &                gather_pattern=p_pat)
    CALL exchange_data(in_array=r_tmp_lat(:,:), out_array=out_lonlat%lat(:), &
      &                gather_pattern=p_pat)

    !-- part 2: exchange vertex lon/lat coordinates:
    DO idim=1,dim3
      CALL exchange_data(in_array=lonv(1:nproma,1:nblks_loc,idim), &
        &                out_array=out_lonlat%lonv(idim,:), &
        &                gather_pattern=p_pat)
      CALL exchange_data(in_array=latv(1:nproma,1:nblks_loc,idim), &
        &                out_array=out_lonlat%latv(idim,:), &
        &                gather_pattern=p_pat)
    END DO ! idim

    ! clean up
    DEALLOCATE(r_tmp_lon, r_tmp_lat, STAT=ierrstat)
    IF (ierrstat /= SUCCESS) CALL finish (routine, 'DEALLOCATE failed.')

  END SUBROUTINE collect_grid_info


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "clonv", "clatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_cells(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, iidx, iblk,                  &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    INTEGER :: max_cell_connectivity
    REAL(wp) :: lonv_temp, latv_temp

    rl_start   = 1
    rl_end     = min_rlcell_int
    i_startblk = p_patch%cells%start_blk(rl_start,1)
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_endblk   = p_patch%cells%end_blk(rl_end,i_nchdom)
    max_cell_connectivity   = p_patch%cells%max_connectivity

!$omp do
    DO jb = i_startblk, i_endblk
      CALL get_indices_c(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO j = 1, max_cell_connectivity
        DO jc = 1, i_startidx - 1
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
        DO jc = i_startidx, i_endidx
          iidx = p_patch%cells%vertex_idx(jc,jb,j)
          iblk = p_patch%cells%vertex_blk(jc,jb,j)
          IF (iidx > 0) THEN
            lonv_temp = p_patch%verts%vertex(iidx,iblk)%lon
            latv_temp = p_patch%verts%vertex(iidx,iblk)%lat
            latv(jc,jb,j) = MERGE(latv_temp, 0.0_wp, &
                 ABS(latv_temp) >= EPSILON(0.0_wp))
            lonv(jc,jb,j) = MERGE(lonv_temp, 0.0_wp, &
                 ABS(lonv_temp) >= EPSILON(0.0_wp))
            IF (ABS(latv_temp) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
              lonv(jc,jb,j) = p_patch%cells%center(jc,jb)%lon
            ENDIF
          ELSE
            lonv(jc,jb,j) = 0.0_wp
            latv(jc,jb,j) = 0.0_wp
          ENDIF
        END DO
        DO jc = i_endidx+1, nproma
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE cf_1_1_grid_cells


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "elonv", "elatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_edges(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER  :: jc, jb, j, iidx, iblk,                 &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    REAL(wp) :: swap(4)
    REAL(wp) :: lonv_temp, latv_temp

    rl_start   = 1
    rl_end     = min_rledge_int
    i_startblk = p_patch%edges%start_blk(rl_start,1)
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_endblk   = p_patch%edges%end_blk(rl_end,i_nchdom)

!$omp do
    DO jb = i_startblk, i_endblk
      CALL get_indices_e(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        iidx = p_patch%edges%vertex_idx(jc,jb,1)
        iblk = p_patch%edges%vertex_blk(jc,jb,1)
        lonv(jc,jb,1) = p_patch%verts%vertex(iidx,iblk)%lon
        latv(jc,jb,1) = p_patch%verts%vertex(iidx,iblk)%lat

        iidx = p_patch%edges%vertex_idx(jc,jb,2)
        iblk = p_patch%edges%vertex_blk(jc,jb,2)
        lonv(jc,jb,3) = p_patch%verts%vertex(iidx,iblk)%lon
        latv(jc,jb,3) = p_patch%verts%vertex(iidx,iblk)%lat

        IF (p_patch%edges%cell_idx(jc,jb,1) > 0) THEN
          iidx = p_patch%edges%cell_idx(jc,jb,1)
          iblk = p_patch%edges%cell_blk(jc,jb,1)
          lonv(jc,jb,4) = p_patch%cells%center(iidx,iblk)%lon
          latv(jc,jb,4) = p_patch%cells%center(iidx,iblk)%lat
        ELSE
          lonv(jc,jb,4) = 0._wp
          latv(jc,jb,4) = 0._wp
        ENDIF

        IF (p_patch%edges%cell_idx(jc,jb,2) > 0) THEN
          iidx = p_patch%edges%cell_idx(jc,jb,2)
          iblk = p_patch%edges%cell_blk(jc,jb,2)
          lonv(jc,jb,2) = p_patch%cells%center(iidx,iblk)%lon
          latv(jc,jb,2) = p_patch%cells%center(iidx,iblk)%lat
        ELSE
          lonv(jc,jb,2) = 0._wp
          latv(jc,jb,2) = 0._wp
        END IF
      END DO
      DO j = 1, 4
        DO jc = 1, i_startidx - 1
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
        DO jc = i_startidx, i_endidx
          latv_temp = MERGE(latv(jc,jb,j), 0.0_wp, &
            &               ABS(latv(jc,jb,j)) >= EPSILON(0.0_wp))
          latv(jc,jb,j) = latv_temp
          IF ( ABS(latv_temp) > 0.5_wp*pi-EPSILON(0.0_wp)) THEN
            lonv(jc,jb,j) = p_patch%edges%center(jc,jb)%lon
          ELSE
            lonv(jc,jb,j) = MERGE(lonv(jc,jb,j), 0.0_wp, &
              &                   ABS(lonv(jc,jb,j)) < EPSILON(0.0_wp))
          END IF
        END DO
        DO jc = i_endidx+1, nproma
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
      END DO
      DO jc = i_startidx, i_endidx
        IF ( check_orientation(p_patch%edges%center(jc,jb)%lon, &
          &                    lonv(jc,jb,:), latv(jc,jb,:), 4) < 0 ) THEN
          swap(1:4) = lonv(jc,jb,4:1:-1)
          lonv(jc,jb,:) = swap(:)
          swap(1:4) = latv(jc,jb,4:1:-1)
          latv(jc,jb,:) = swap(:)
        END IF
      END DO
    END DO
!$omp end do nowait

  END SUBROUTINE cf_1_1_grid_edges


  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "vlonv", "vlatv" into a CF1.1 compatible form for ocean
  ! This is needed for general connectivity grids, the original routine  cf_1_1_grid_verts will not work
  !
  ! This is also sloppy for boundary vertices, as they seem not to be defined
  !
  !  based on SUBROUTINE cf_1_1_grid_verts
  !
  SUBROUTINE cf_1_1_grid_verts_ocean(patch_2D, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: patch_2D
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, &
      &        iidx(nproma,patch_2D%verts%max_connectivity), &
      &        iblk(nproma,patch_2D%verts%max_connectivity), &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    INTEGER :: last_valid_cell(nproma), max_vrtx_conn

    rl_start   = 2
    rl_end     = min_rlvert
    i_nchdom   = MAX(1,patch_2D%n_childdom)
    i_startblk = patch_2D%verts%start_blk(rl_start,1)
    i_endblk   = patch_2D%verts%end_blk(rl_end,i_nchdom)
    max_vrtx_conn = patch_2D%verts%max_connectivity

!$omp do
    DO jb = i_startblk, i_endblk
      CALL get_indices_v(patch_2D, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      last_valid_cell = 0
      DO j = 1, max_vrtx_conn
        DO jc = 1, i_startidx - 1
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
        DO jc = i_startidx, i_endidx
          last_valid_cell(jc) = MERGE(j, last_valid_cell(jc), &
            &                         patch_2D%verts%cell_idx(jc,jb,j) > 0)
          iidx(jc,j) = patch_2D%verts%cell_idx(jc,jb,last_valid_cell(jc))
          iblk(jc,j) = patch_2D%verts%cell_blk(jc,jb,last_valid_cell(jc))
        END DO
        DO jc = i_startidx, i_endidx
          lonv(jc,jb, j) = patch_2D%cells%center(iidx(jc,j),iblk(jc,j))%lon
          latv(jc,jb, j) = patch_2D%cells%center(iidx(jc,j),iblk(jc,j))%lat
        ENDDO
        DO jc = i_endidx+1, nproma
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
      END DO
    END DO
!$omp end do nowait
  END SUBROUTINE cf_1_1_grid_verts_ocean
  !------------------------------------------------------------------------------------------------



  !------------------------------------------------------------------------------------------------
  !> Reformat vertex lon/lat coordinates "vlonv", "vlatv" into a CF1.1 compatible form.
  !
  !  based on SUBROUTINE mo_gridrefinement::write_patch(p)
  !
  SUBROUTINE cf_1_1_grid_verts(p_patch, lonv, latv)
    TYPE(t_patch),      INTENT(IN)    :: p_patch
    REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
    ! local variables
    INTEGER :: jc, jb, j, iidx, iblk,                  &
      &        rl_start, rl_end, i_startblk, i_endblk, &
      &        i_startidx, i_endidx, i_nchdom
    INTEGER :: max_vertex_connectivity

    rl_start   = 2
    rl_end     = min_rlvert
    i_nchdom   = MAX(1,p_patch%n_childdom)
    i_startblk = p_patch%verts%start_blk(rl_start,1)
    i_endblk   = p_patch%verts%end_blk(rl_end,i_nchdom)
    max_vertex_connectivity = p_patch%verts%max_connectivity

!$omp do
    DO jb = i_startblk, i_endblk
      CALL get_indices_v(p_patch, jb, i_startblk, i_endblk, &
        &                i_startidx, i_endidx, rl_start, rl_end)
      DO j = 1,max_vertex_connectivity
        DO jc = 1, i_startidx - 1
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
        DO jc = i_startidx, i_endidx
          IF ((p_patch%verts%cell_idx(jc,jb,j) == 0) .AND. &
            & (p_patch%verts%refin_ctrl(jc,jb) /= 1)) THEN
            iidx = p_patch%verts%cell_idx(jc,jb,5)
            iblk = p_patch%verts%cell_blk(jc,jb,5)
            lonv(jc,jb,max_vertex_connectivity+1-j) = p_patch%cells%center(iidx,iblk)%lon
            latv(jc,jb,max_vertex_connectivity+1-j) = p_patch%cells%center(iidx,iblk)%lat
          ELSE IF ((p_patch%verts%cell_idx(jc,jb,j) < 0) .OR. &
            &      (p_patch%verts%refin_ctrl(jc,jb) == 1)) THEN
            lonv(jc,jb,max_vertex_connectivity+1-j) = 0._wp
            latv(jc,jb,max_vertex_connectivity+1-j) = 0._wp
          ELSE
            iidx = p_patch%verts%cell_idx(jc,jb,j)
            iblk = p_patch%verts%cell_blk(jc,jb,j)
            lonv(jc,jb,max_vertex_connectivity+1-j) = p_patch%cells%center(iidx,iblk)%lon
            latv(jc,jb,max_vertex_connectivity+1-j) = p_patch%cells%center(iidx,iblk)%lat
          ENDIF
        ENDDO
        DO jc = i_endidx+1, nproma
          lonv(jc,jb,j) = 0.0_wp
          latv(jc,jb,j) = 0.0_wp
        END DO
      ENDDO
    END DO
!$omp end do nowait
  END SUBROUTINE cf_1_1_grid_verts


  !------------------------------------------------------------------------------------------------
  !> Sets the grid information in output file
  !
  SUBROUTINE set_grid_info_netcdf(gridID, grid_info)
    INTEGER,             INTENT(IN) :: gridID
    TYPE(t_grid_info),   INTENT(IN) :: grid_info
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_grid_info_netcdf"

    CALL gridDefXvals  (gridID, grid_info%lon)
    CALL gridDefYvals  (gridID, grid_info%lat)
    CALL gridDefXbounds(gridID, grid_info%lonv)
    CALL gridDefYbounds(gridID, grid_info%latv)
  END SUBROUTINE set_grid_info_netcdf


  !------------------------------------------------------------------------------------------------
  !> Depending on the vertical interpolation type of an output file,
  !  this function returns the list of variable names (ml_varlist,
  !  pl_varlist, etc.)
  !
  SUBROUTINE get_varlist(of, p_varlist)
    TYPE (t_output_file), INTENT(INOUT) :: of
    CHARACTER(LEN=vname_len), POINTER :: p_varlist(:)

    SELECT CASE(of%ilev_type)
    CASE(level_type_ml)
      p_varlist => of%name_list%ml_varlist
    CASE(level_type_pl)
      p_varlist => of%name_list%pl_varlist
    CASE(level_type_hl)
      p_varlist => of%name_list%hl_varlist
    CASE(level_type_il)
      p_varlist => of%name_list%il_varlist
    END SELECT
  END SUBROUTINE get_varlist


  !------------------------------------------------------------------------------------------------
  !> Declaration of the grid information (RLAT/RLON) in output file, GRIB2 format.
  !
  SUBROUTINE set_grid_info_grb2(of)
    TYPE (t_output_file), INTENT(INOUT) :: of
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::set_grid_info_grb"
    CHARACTER(LEN=4), PARAMETER :: grid_coord_name(2) = (/ "RLON", "RLAT" /)
    TYPE (t_grib2_var) :: grid_coord_grib2(2)
    INTEGER :: igrid,i,vlistID,gridID(3),zaxisID
    TYPE(t_verticalAxis), POINTER  :: zaxis
    INTEGER, PARAMETER :: idx(3) = (/ ICELL, IEDGE, IVERT /)

    ! geographical longitude RLON
    grid_coord_grib2(1) = grib2_var(               0,   &  ! discipline
      &                                          191,   &  ! category
      &                                            2,   &  ! number
      &                              DATATYPE_PACK16,   &  ! bits
      &                            GRID_UNSTRUCTURED,   &  ! gridtype
      &                                    GRID_CELL )     ! subgridtype

    ! geographical latitude RLAT
    grid_coord_grib2(2) = grib2_var(               0,   &  ! discipline
      &                                          191,   &  ! category
      &                                            1,   &  ! number
      &                              DATATYPE_PACK16,   &  ! bits
      &                            GRID_UNSTRUCTURED,   &  ! gridtype
      &                                    GRID_CELL )     ! subgridtype

    vlistID =  of%cdiVlistID
    zaxis => of%verticalAxisList%getEntry(icon_zaxis_type=ZA_surface)
    IF (.NOT. ASSOCIATED(zaxis))  CALL finish(routine, 'Zaxis undefined.')
    zaxisID = zaxis%cdi_id

    SELECT CASE(of%name_list%remap)
    CASE (REMAP_NONE)
      gridID(1) = of%cdiCellGridID
      gridID(2) = of%cdiEdgeGridID
      gridID(3) = of%cdiVertGridID
      DO igrid=1,3
        DO i=1,2 ! for longitude, latitude:
          of%cdi_grb2(idx(igrid),i) = vlistDefVar(vlistID, gridID(igrid), zaxisID, TIME_CONSTANT)
          CALL vlistDefVarDatatype(vlistID,  of%cdi_grb2(idx(igrid),i), grid_coord_grib2(i)%bits)
          CALL vlistDefVarTsteptype(vlistID, of%cdi_grb2(idx(igrid),i), TSTEP_CONSTANT)
          CALL vlistDefVarName(vlistID,      of%cdi_grb2(idx(igrid),i), TRIM(grid_coord_name(i)))

          ! Set GRIB2 Triplet
          CALL vlistDefVarParam( vlistID, of%cdi_grb2(idx(igrid),i),  &
            &  cdiEncodeParam(grid_coord_grib2(i)%number,             &
            &                 grid_coord_grib2(i)%category,           &
            &                 grid_coord_grib2(i)%discipline) )

          ! GRIB2 Quick hack: Set additional GRIB2 keys
          CALL set_GRIB2_additional_keys(vlistID, of%cdi_grb2(idx(igrid),i),   &
            &                            gribout_config(of%phys_patch_id) )

          ! Set ensemble keys in SECTION 4 (if applicable)
          CALL set_GRIB2_ensemble_keys(vlistID, of%cdi_grb2(idx(igrid),i),   &
            &                            gribout_config(of%phys_patch_id) )

          ! Set local use SECTION 2
          CALL set_GRIB2_local_keys(vlistID, of%cdi_grb2(idx(igrid),i),   &
            &                       gribout_config(of%phys_patch_id) )
        END DO
      END DO

    CASE (REMAP_REGULAR_LATLON)
      ! for longitude, latitude:
      DO i=1,2
        of%cdi_grb2(ILATLON,i) = vlistDefVar(vlistID, of%cdiLonLatGridID, zaxisID, TIME_CONSTANT)
        CALL vlistDefVarDatatype(vlistID,  of%cdi_grb2(ILATLON,i), grid_coord_grib2(i)%bits)
        CALL vlistDefVarTsteptype(vlistID, of%cdi_grb2(ILATLON,i), TSTEP_CONSTANT)
        CALL vlistDefVarName(vlistID,      of%cdi_grb2(ILATLON,i), TRIM(grid_coord_name(i)))
        ! Set GRIB2 Triplet
        CALL vlistDefVarParam( vlistID, of%cdi_grb2(ILATLON,i),   &
          &  cdiEncodeParam(grid_coord_grib2(i)%number,            &
          &                 grid_coord_grib2(i)%category,          &
          &                 grid_coord_grib2(i)%discipline) )

        ! GRIB2 Quick hack: Set additional GRIB2 keys
        CALL set_GRIB2_additional_keys(vlistID, of%cdi_grb2(ILATLON,i), &
          &                            gribout_config(of%phys_patch_id) )

        ! Set ensemble keys in SECTION 4 (if applicable)
        CALL set_GRIB2_ensemble_keys(vlistID, of%cdi_grb2(ILATLON,i), &
          &                            gribout_config(of%phys_patch_id) )

        ! Set local use SECTION 2
        CALL set_GRIB2_local_keys(vlistID, of%cdi_grb2(ILATLON,i),    &
          &                       gribout_config(of%phys_patch_id) )
      END DO

    CASE DEFAULT
      CALL finish(routine, "Unsupported grid type.")
    END SELECT

  END SUBROUTINE set_grid_info_grb2


  !------------------------------------------------------------------------------------------------
  !
  ! Copies the grid information from grid file to output file
  !
  SUBROUTINE copy_grid_info(of, patch_info)
    TYPE(t_patch_info),   INTENT(INOUT) :: patch_info
    TYPE (t_output_file), INTENT(IN) :: of

    INTEGER :: ncid, dimid, varid, tlen
    INTEGER :: i_nc, i_ne, i_nv, max_cell_connectivity, max_verts_connectivity

    REAL(wp), ALLOCATABLE :: clon(:), clat(:), clonv(:,:), clatv(:,:)
    REAL(wp), ALLOCATABLE :: elon(:), elat(:), elonv(:,:), elatv(:,:)
    REAL(wp), ALLOCATABLE :: vlon(:), vlat(:), vlonv(:,:), vlatv(:,:)

    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::copy_grid_info"

    ! Please note: The following is more or less a copy from mo_io_vlist with adaptions
    ! to the data structures used here.
    ! Unfortunately it seems necessary to open the gridfile for reading the information
    ! since it is not read and stored during patch input.

    !---------------------------------------------------------------------------
    ! Open grid file, read dimensions and make a cross check if they match.
    ! This is just for safety and could be skipped, of course.

    tlen = LEN_TRIM(patch_info%grid_filename)
    CALL nf(nf_open(patch_info%grid_filename(1:tlen), NF_NOWRITE, ncid))

    CALL nf(nf_inq_dimid(ncid, 'nv', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, max_cell_connectivity))
    CALL nf(nf_inq_dimid(ncid, 'ne', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, max_verts_connectivity))
    !
    CALL nf(nf_inq_dimid(ncid, 'cell', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_nc))
    !
    CALL nf(nf_inq_dimid(ncid, 'edge', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_ne))
    !
    CALL nf(nf_inq_dimid(ncid, 'vertex', dimid))
    CALL nf(nf_inq_dimlen(ncid, dimid, i_nv))

    IF(i_nc /= patch_info%grid_info(icell)%n_log) &
      CALL finish(routine,'Number of cells differs in '//patch_info%grid_filename(1:tlen))
    IF(i_ne /= patch_info%grid_info(iedge)%n_log) &
      CALL finish(routine,'Number of edges differs in '//patch_info%grid_filename(1:tlen))
    IF(i_nv /= patch_info%grid_info(ivert)%n_log) &
      CALL finish(routine,'Number of verts differs in '//patch_info%grid_filename(1:tlen))
    !
    !---------------------------------------------------------------------------
    ! cell grid
    CALL nf(nf_inq_varid(ncid, 'clon', varid))
    ALLOCATE(clon(i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clon))
    CALL reorder1(patch_info%grid_info(icell)%log_dom_starts, &
      &           patch_info%grid_info(icell)%log_dom_counts,clon)
    CALL gridDefXvals(of%cdiCellGridID, clon)
    DEALLOCATE(clon)


    CALL nf(nf_inq_varid(ncid, 'clat', varid))
    ALLOCATE(clat(i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clat))
    CALL reorder1(patch_info%grid_info(icell)%log_dom_starts, &
      &           patch_info%grid_info(icell)%log_dom_counts,clat)

    CALL gridDefYvals(of%cdiCellGridID, clat)
    DEALLOCATE(clat)


    CALL nf(nf_inq_varid(ncid, 'clon_vertices', varid))
    ALLOCATE(clonv(max_cell_connectivity, i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clonv))
    CALL reorder2(patch_info%grid_info(icell)%log_dom_starts, &
      &           patch_info%grid_info(icell)%log_dom_counts,clonv)

    CALL gridDefXbounds(of%cdiCellGridID, clonv)
    DEALLOCATE(clonv)


    CALL nf(nf_inq_varid(ncid, 'clat_vertices', varid))
    ALLOCATE(clatv(max_cell_connectivity, i_nc))
    CALL nf(nf_get_var_double(ncid, varid, clatv))
    CALL reorder2(patch_info%grid_info(icell)%log_dom_starts, &
      &           patch_info%grid_info(icell)%log_dom_counts,clatv)

    CALL gridDefYbounds(of%cdiCellGridID, clatv)
    DEALLOCATE(clatv)

    !-------------------------------------------------------------------------
    ! edge grid

    ALLOCATE(elon(i_ne))
    CALL nf(nf_inq_varid(ncid, 'elon', varid))
    CALL nf(nf_get_var_double(ncid, varid, elon))
    CALL reorder1(patch_info%grid_info(iedge)%log_dom_starts, &
      &           patch_info%grid_info(iedge)%log_dom_counts,elon)

    CALL gridDefXvals(of%cdiEdgeGridID, elon)
    DEALLOCATE(elon)

    ALLOCATE(elat(i_ne))
    CALL nf(nf_inq_varid(ncid, 'elat', varid))
    CALL nf(nf_get_var_double(ncid, varid, elat))
    CALL reorder1(patch_info%grid_info(iedge)%log_dom_starts, &
      &           patch_info%grid_info(iedge)%log_dom_counts,elat)

    CALL gridDefYvals(of%cdiEdgeGridID, elat)
    DEALLOCATE(elat)

    ALLOCATE(elonv(4, i_ne))
    CALL nf(nf_inq_varid(ncid, 'elon_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elonv))
    CALL reorder2(patch_info%grid_info(iedge)%log_dom_starts, &
      &           patch_info%grid_info(iedge)%log_dom_counts,elonv)

    CALL gridDefXbounds(of%cdiEdgeGridID, elonv)
    DEALLOCATE(elonv)

    ALLOCATE(elatv(4, i_ne))
    CALL nf(nf_inq_varid(ncid, 'elat_vertices', varid))
    CALL nf(nf_get_var_double(ncid, varid, elatv))
    CALL reorder2(patch_info%grid_info(iedge)%log_dom_starts, &
      &           patch_info%grid_info(iedge)%log_dom_counts,elatv)

    CALL gridDefYbounds(of%cdiEdgeGridID, elatv)
    DEALLOCATE(elatv)

    !-------------------------------------------------------------------------
    ! vertex grid
    CALL nf(nf_inq_varid(ncid, 'vlon', varid))
    ALLOCATE(vlon(i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlon))
    CALL reorder1(patch_info%grid_info(ivert)%log_dom_starts, &
      &           patch_info%grid_info(ivert)%log_dom_counts,vlon)

    CALL gridDefXvals(of%cdiVertGridID, vlon)
    DEALLOCATE(vlon)

    CALL nf(nf_inq_varid(ncid, 'vlat', varid))
    ALLOCATE(vlat(i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlat))
    CALL reorder1(patch_info%grid_info(ivert)%log_dom_starts, &
      &           patch_info%grid_info(ivert)%log_dom_counts,vlat)

    CALL gridDefYvals(of%cdiVertGridID, vlat)
    DEALLOCATE(vlat)

    CALL nf(nf_inq_varid(ncid, 'vlon_vertices', varid))
    ALLOCATE(vlonv(max_verts_connectivity, i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlonv))
    CALL reorder2(patch_info%grid_info(ivert)%log_dom_starts, &
      &           patch_info%grid_info(ivert)%log_dom_counts,vlonv)

    CALL gridDefXbounds(of%cdiVertGridID, vlonv)
    DEALLOCATE(vlonv)

    CALL nf(nf_inq_varid(ncid, 'vlat_vertices', varid))
    ALLOCATE(vlatv(max_verts_connectivity, i_nv))
    CALL nf(nf_get_var_double(ncid, varid, vlatv))
    CALL reorder2(patch_info%grid_info(ivert)%log_dom_starts, &
      &           patch_info%grid_info(ivert)%log_dom_counts, vlatv)

    CALL gridDefYbounds(of%cdiVertGridID, vlatv)
    DEALLOCATE(vlatv)

    !-------------------------------------------------------------------------

    ! Close NetCDF file, it is not needed any more
    CALL nf(nf_close(ncid))

  CONTAINS

    SUBROUTINE nf(errstat)
      USE mo_netcdf_errhandler, ONLY: nf_ext => nf
      INTEGER, INTENT(IN) :: errstat

      CALL nf_ext(errstat, routine)
    END SUBROUTINE nf

    ! reorder1: get the physical patch points from the logical patch
    ! Note that this works within the array as long as idx is monotonically increasing
    SUBROUTINE reorder1(starts, counts, array)
      INTEGER, INTENT(IN), CONTIGUOUS :: starts(:), counts(:)
      REAL(wp), INTENT(INOUT), CONTIGUOUS :: array(:)
      INTEGER :: i, j, m, n, doff, soff

      m = SIZE(counts)
      doff = 0
      DO i = 1, m
        n = counts(i)
        soff = starts(i) - 1
        DO j = 1, n
          array(doff+j) = array(soff+j)
        END DO
        doff = doff + n
      ENDDO
    END SUBROUTINE reorder1

    ! reorder2: same as reorder1 for 2D array
    SUBROUTINE reorder2(starts, counts, array)
      INTEGER, INTENT(IN), CONTIGUOUS :: starts(:), counts(:)
      REAL(wp), INTENT(INOUT), CONTIGUOUS :: array(:,:)
      INTEGER :: i, j, m, n, doff, soff

      m = SIZE(counts)
      doff = 0
      DO i = 1, m
        n = counts(i)
        soff = starts(i) - 1
        DO j = 1, n
          array(:,doff+j) = array(:,soff+j)
        END DO
        doff = doff + n
      ENDDO
    END SUBROUTINE reorder2

  END SUBROUTINE copy_grid_info

  SUBROUTINE allgather_grid_info_cve(connectivity, nblks, nblks_glb, &
    &                                coordinates, grid_info, gather_pattern, &
    &                                cf_1_1_grid, keep_grid_info, p_patch)
    INTEGER,            INTENT(IN)                     :: connectivity
    INTEGER,            INTENT(IN)                     :: nblks
    INTEGER,            INTENT(IN)                     :: nblks_glb
    TYPE(t_geographical_coordinates), INTENT(IN)       :: coordinates(:,:)
    TYPE(t_grid_info), TARGET, INTENT(INOUT)           :: grid_info
    TYPE(t_comm_gather_pattern), TARGET, INTENT(INOUT) :: gather_pattern
    INTERFACE
      SUBROUTINE cf_1_1_grid(p_patch, lonv, latv)
        USE mo_model_domain, ONLY: t_patch
        USE mo_kind,         ONLY: wp
        TYPE(t_patch),      INTENT(IN)    :: p_patch
        REAL(wp),           INTENT(INOUT) :: lonv(:,:,:), latv(:,:,:)
      END SUBROUTINE cf_1_1_grid
    END INTERFACE
    !> only those I/O processes which need the coordinate data for
    !! this patch set keep_grid_info to .TRUE.
    LOGICAL,            INTENT(IN)                     :: keep_grid_info
    TYPE(t_patch), OPTIONAL, INTENT(IN)                :: p_patch

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allgather_grid_info_cve"
    INTEGER :: i
    LOGICAL :: is_io
    TYPE(t_comm_allgather_pattern) :: allgather_pattern
    REAL(wp), ALLOCATABLE          :: lonv(:,:,:), latv(:,:,:)
    REAL(wp), POINTER              :: r1d(:)
    REAL(wp), TARGET               :: dummy(1)

    is_io = my_process_is_io()
    IF (is_io .AND. keep_grid_info) THEN
      ALLOCATE(grid_info%lon(nproma*nblks_glb),                &
        &      grid_info%lat(nproma*nblks_glb),                &
        &      grid_info%lonv(connectivity, nproma*nblks_glb), &
        &      grid_info%latv(connectivity, nproma*nblks_glb), &
        &      lonv(1,1,connectivity), latv(1,1,connectivity))
      grid_info%lon = 0._wp
      grid_info%lat = 0._wp
      grid_info%lonv = 0._wp
      grid_info%latv = 0._wp
    ELSE IF (is_io .AND. .NOT. keep_grid_info) THEN
      ALLOCATE(r1d(nproma*nblks_glb), lonv(0,0,connectivity), latv(0,0,connectivity))
    ELSE
      ALLOCATE(lonv(nproma, nblks, connectivity), &
        &      latv(nproma, nblks, connectivity))
!$omp parallel
      CALL cf_1_1_grid(p_patch, lonv, latv)
!$omp end parallel
      r1d => dummy
    END IF

    CALL setup_comm_allgather_pattern(gather_pattern, p_comm_work_2_io, &
      &                               allgather_pattern)

    ! gathers coordinates on all io procs
    IF (is_io .AND. keep_grid_info) r1d => grid_info%lon
    CALL exchange_data(in_array=coordinates(:,:)%lon, out_array=r1d, &
      &                allgather_pattern=allgather_pattern)
    IF (is_io .AND. keep_grid_info) r1d => grid_info%lat
    CALL exchange_data(in_array=coordinates(:,:)%lat, out_array=r1d, &
      &                allgather_pattern=allgather_pattern)
    DO i = 1, connectivity
      IF (is_io .AND. keep_grid_info) r1d => grid_info%lonv(i,:)
      CALL exchange_data(in_array=lonv(:,:,i), out_array=r1d, &
        &                allgather_pattern=allgather_pattern)
      IF (is_io .AND. keep_grid_info) r1d => grid_info%latv(i,:)
      CALL exchange_data(in_array=latv(:,:,i), out_array=r1d, &
        &                allgather_pattern=allgather_pattern)
    END DO

    IF (is_io .AND. .NOT. keep_grid_info) DEALLOCATE(r1d)

    CALL delete_comm_allgather_pattern(allgather_pattern)

    DEALLOCATE (lonv, latv)
  END SUBROUTINE allgather_grid_info_cve

  SUBROUTINE allgather_grid_info(patch_info, keep_grid_info, p_patch)
    TYPE(t_patch_info),    INTENT(INOUT) :: patch_info
    !> only those I/O processes which need the coordinate data for
    !! this patch set keep_grid_info to .TRUE.
    LOGICAL,               INTENT(IN)    :: keep_grid_info
    TYPE(t_patch), OPTIONAL, TARGET, INTENT(IN)    :: p_patch

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::allgather_grid_info"
    TYPE(t_geographical_coordinates), TARGET :: dummy_coordinates(1,1)
    TYPE(t_geographical_coordinates), POINTER :: coordinates(:,:)
    INTEGER :: dummy_global_size, dummy_owner_local(0), dummy_glb_index(0)
    TYPE(t_comm_gather_pattern), TARGET :: empty_gather_pattern
    TYPE(t_comm_gather_pattern), POINTER :: gather_pattern
    LOGICAL :: is_io
    INTEGER :: nblks

    is_io = my_process_is_io()
    IF (is_io) THEN
      dummy_global_size = 0
      CALL setup_comm_gather_pattern(dummy_global_size, dummy_owner_local, &
        &                            dummy_glb_index, empty_gather_pattern)
      dummy_coordinates(1,1)%lon = -1._wp
      dummy_coordinates(1,1)%lat = -1._wp
      coordinates => dummy_coordinates
      gather_pattern => empty_gather_pattern
      nblks = 0
    ELSE
      coordinates => p_patch%cells%center
      gather_pattern => patch_info%p_pat_c
      nblks = p_patch%nblks_c
    END IF

    CALL allgather_grid_info_cve(patch_info%max_cell_connectivity, &
      &                          nblks, patch_info%nblks_glb_c, &
      &                          coordinates, patch_info%grid_info(icell), &
      &                          gather_pattern, cf_1_1_grid_cells, &
      &                          keep_grid_info, p_patch)

    IF (.NOT. is_io) THEN
      coordinates => p_patch%verts%vertex
      gather_pattern => patch_info%p_pat_v
      nblks = p_patch%nblks_v
    END IF

    IF (my_process_is_oceanic()) THEN
      CALL allgather_grid_info_cve(patch_info%max_vertex_connectivity, &
        &                          nblks, patch_info%nblks_glb_v, &
        &                          coordinates, patch_info%grid_info(ivert), &
        &                          gather_pattern, cf_1_1_grid_verts_ocean, &
        &                          keep_grid_info, p_patch)
    ELSE
      CALL allgather_grid_info_cve(9 - patch_info%max_cell_connectivity, &
        &                          nblks, patch_info%nblks_glb_v, &
        &                          coordinates, patch_info%grid_info(ivert), &
        &                          gather_pattern, cf_1_1_grid_verts, &
        &                          keep_grid_info, p_patch)
    END IF

    IF (.NOT. is_io) THEN
      coordinates => p_patch%edges%center
      gather_pattern => patch_info%p_pat_e
      nblks = p_patch%nblks_e
    END IF

    CALL allgather_grid_info_cve(4, &
      &                          nblks, patch_info%nblks_glb_e, &
      &                          coordinates, patch_info%grid_info(iedge), &
      &                          gather_pattern, cf_1_1_grid_edges, &
      &                          keep_grid_info, p_patch)

    IF(is_io) CALL delete_comm_gather_pattern(empty_gather_pattern)

  END SUBROUTINE allgather_grid_info


  !------------------------------------------------------------------------------------------------
  !> Writes the grid information in output file, GRIB2 format.
  !
! #ifndef __NO_ICON_ATMO__
  SUBROUTINE write_grid_info_grb2(of, patch_info)
    TYPE (t_output_file), INTENT(INOUT)           :: of
    TYPE(t_patch_info),   INTENT(IN),   TARGET    :: patch_info (:)
    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//"::write_grid_info_grb2"

    INTEGER                        :: errstat, idom, igrid, n, idom_log
    TYPE (t_lon_lat_grid), POINTER :: grid
    REAL(wp), ALLOCATABLE          :: rotated_pts(:,:,:), r_out_dp(:,:), r_out_dp_1D(:)
    CHARACTER(LEN=vname_len), POINTER :: p_varlist(:)
    INTEGER, PARAMETER :: idx(3) = (/ ICELL, IEDGE, IVERT /)

    CALL get_varlist(of, p_varlist)
    SELECT CASE(of%name_list%remap)
    CASE (REMAP_NONE)
      idom     = of%phys_patch_id
      idom_log = patch_info(idom)%log_patch_id
      DO igrid=1,3
#ifdef HAVE_CDI_PIO
        IF (pio_type == pio_type_cdipio) THEN
          n = SIZE(patch_info(idom_log)%grid_info(igrid)%lon)
        ELSE
#endif
          n = patch_info(idom_log)%ri(igrid)%n_glb
#ifdef HAVE_CDI_PIO
        END IF
#endif
        ! allocate data buffer:
        ALLOCATE(r_out_dp_1D(n), stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed!')

        ! write RLON, RLAT
        IF (one_of(GRB2_GRID_INFO_NAME(igrid,IRLON), p_varlist) /= -1) THEN
          r_out_dp_1D(:) = patch_info(idom_log)%grid_info(igrid)%lon(1:n) * rad2deg
          CALL write_unstruct_grid2var(of%cdiFileID, &
            of%cdi_grb2(idx(igrid),IRLON), r_out_dp_1D, &
            patch_info(idom_log)%ri(igrid))
        END IF
        IF (one_of(GRB2_GRID_INFO_NAME(igrid,IRLAT), p_varlist) /= -1) THEN
          r_out_dp_1D(:) = patch_info(idom_log)%grid_info(igrid)%lat(1:n) * rad2deg
          CALL write_unstruct_grid2var(of%cdiFileID, &
            of%cdi_grb2(idx(igrid),IRLAT), r_out_dp_1D, &
            patch_info(idom_log)%ri(igrid))
        END IF

        ! clean up
        DEALLOCATE(r_out_dp_1D, stat=errstat)
        IF (errstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed!')
      END DO

    CASE (REMAP_REGULAR_LATLON)
      ! allocate data buffer:
      grid => lonlat_grids%list(of%name_list%lonlat_id)%grid
      ! compute some entries of lon-lat grid specification:
      CALL compute_lonlat_specs(grid)
      ALLOCATE(rotated_pts(grid%lon_dim, grid%lat_dim, 2), &
        &      r_out_dp(grid%lon_dim,grid%lat_dim), stat=errstat)
      IF (errstat /= SUCCESS) CALL finish(routine, 'ALLOCATE failed!')
      ! compute grid points of rotated lon/lat grid
      CALL rotate_latlon_grid(grid, rotated_pts)

      ! write RLON, RLAT
      IF (one_of(GRB2_GRID_INFO_NAME(0,IRLON), p_varlist) /= -1) THEN
        r_out_dp(:,:) = rotated_pts(:,:,1) * rad2deg
        CALL write_remap_grid2var(of%cdiFileID, of%cdi_grb2(ILATLON,IRLON), &
          r_out_dp)
      END IF
      IF (one_of(GRB2_GRID_INFO_NAME(0,IRLAT), p_varlist) /= -1) THEN
        r_out_dp(:,:) = rotated_pts(:,:,2) * rad2deg
        CALL write_remap_grid2var(of%cdiFileID, of%cdi_grb2(ILATLON,IRLAT), &
          r_out_dp)
      END IF

      ! clean up
      DEALLOCATE(rotated_pts, r_out_dp, stat=errstat)
      IF (errstat /= SUCCESS) CALL finish(routine, 'DEALLOCATE failed!')

    CASE DEFAULT
      CALL finish(routine, "Unsupported grid type.")
    END SELECT

  CONTAINS
    SUBROUTINE write_unstruct_grid2var(fileid, varid, r_out_dp_1D, ri)
      INTEGER, INTENT(in) :: fileid, varid
      REAL(wp), INTENT(in), CONTIGUOUS :: r_out_dp_1d(:)
      TYPE(t_reorder_info), INTENT(in) :: ri
#ifdef HAVE_CDI_PIO
      TYPE(extent) :: grid_size_desc, grid_part_desc
      INTEGER(c_int) :: grid_chunk(2, 3)
      INTEGER :: i,j
#endif

#ifdef HAVE_CDI_PIO
      IF (pio_type == pio_type_cdipio) THEN
        CALL streamWriteVarPart(fileid, varid, r_out_dp_1D(i:j), 0, &
          &                     ri%reorder_idxlst_xt(1))
      ELSE
#endif
        CALL streamWriteVar(fileid, varid, r_out_dp_1D, 0)
#ifdef HAVE_CDI_PIO
      END IF
#endif
    END SUBROUTINE write_unstruct_grid2var

    SUBROUTINE write_remap_grid2var(fileid, varid, r_out_dp)
      INTEGER, INTENT(in) :: fileid, varid
      REAL(wp), TARGET, INTENT(in), CONTIGUOUS :: r_out_dp(:,:)
#ifdef HAVE_CDI_PIO
      TYPE(extent) :: grid_size_desc(2), grid_part_desc(2)
      INTEGER(c_int) :: grid_chunk(2, 3)
      INTEGER :: p(2), q(2), div(2), div2, partx, party, i
#endif

#ifdef HAVE_CDI_PIO
      IF (pio_type == pio_type_cdipio) THEN
        div2 = FLOOR(SQRT(REAL(p_n_work)))
        DO WHILE (div2 > 2)
          IF (MOD(p_n_work, div2) == 0) EXIT
          div2 = div2 - 1
        END DO
        div(1) = p_n_work / div2
        div(2) = div2
        p(1) = MOD(p_pe_work, div(1)) + 1
        p(2) = p_pe_work / div(1) + 1
        DO i = 1, 2
          grid_size_desc(i)%first = 1
          grid_size_desc(i)%size = SIZE(r_out_dp, i)
          grid_part_desc(i) = uniform_partition(grid_size_desc(i), div(i), p(i))
          p(i) = grid_part_desc(i)%first
          grid_chunk(i,1) = INT(p(i), c_int)
          q(i) = p(i) + grid_part_desc(i)%size - 1
          grid_chunk(i,2) = INT(q(i), c_int)
        END DO
        grid_chunk(1,3) = 1
        grid_chunk(2,3) = 1
        CALL streamWriteVarChunk(fileid, varid, grid_chunk, &
          &                      r_out_dp(p(1):q(1),p(2):q(2)), 0)
      ELSE
#endif
        CALL streamWriteVar(fileid, varid, r_out_dp, 0)
#ifdef HAVE_CDI_PIO
      END IF
#endif
    END SUBROUTINE write_remap_grid2var

  END SUBROUTINE write_grid_info_grb2
! #endif

END MODULE mo_name_list_output_gridinfo
