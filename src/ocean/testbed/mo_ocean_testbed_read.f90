! @brief tests the mo_read_interface methods
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

MODULE mo_ocean_testbed_read

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, message_text, finish
  USE mo_mpi,                 ONLY: work_mpi_barrier, my_process_is_stdio, &
    &                               p_n_work, p_pe_work, p_comm_work
  USE mo_timer,               ONLY: init_timer, ltimer, new_timer, timer_start, timer_stop, &
    & activate_sync_timers, timers_level, timer_barrier, timer_radiaton_recv
  USE mo_parallel_config,     ONLY: nproma, icon_comm_method
  USE mo_ocean_nml
  USE mo_master_control,      ONLY: get_my_process_name, get_my_model_no
  USE mo_io_units,            ONLY: filename_max
  USE mo_io_config,           ONLY:  read_netcdf_broadcast_method, &
    &                                read_netcdf_distribute_method

  USE mo_model_domain,        ONLY: t_patch, t_patch_3D
  USE mo_grid_config,         ONLY: dynamics_grid_filename
!   USE mo_test_netcdf_read,    ONLY: netcdf_write_oncells_3D_time
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf
  USE mo_read_interface, ONLY: t_stream_id, read_2D, on_cells, on_vertices, on_edges, &
    & read_2d_int, read_3d, read_2d_time, read_2d_extdim, read_2d_extdim_int, &
    & read_3d_time, read_3d_extdim, read_2d_1time, read_2d_1lev_1time, read_3d_1time, &
    & openinputfile, closefile
  USE mo_read_netcdf_distributed, ONLY: setup_distrib_read, delete_distrib_read
  USE mo_communication,       ONLY: idx_no, blk_no, makeScatterPattern
  USE mo_util_file,           ONLY: util_unlink

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

  PUBLIC :: ocean_test_read

CONTAINS

  !-------------------------------------------------------------------------
  !>
  !! Test reading T
  SUBROUTINE ocean_test_read(namelist_filename, shr_namelist_filename, patch_3d)

    USE iso_c_binding, ONLY : c_int

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename
    TYPE(t_patch_3d ),TARGET, INTENT(in)    :: patch_3d

    REAL(wp), POINTER :: T(:,:,:,:), T_check(:,:,:,:)   ! is (nproma, levels, blocks, time )
    INTEGER :: levels, lnwl_size, return_status
    INTEGER(c_int) :: return_status_c
    REAL(wp), POINTER :: cell_data_broadcast(:,:), cell_data_distribute(:,:) ! is (nproma, blocks )
    REAL(wp), POINTER :: vertex_data_broadcast(:,:), vertex_data_distribute(:,:) ! is (nproma, blocks )
    REAL(wp), POINTER :: edge_data_broadcast(:,:), edge_data_distribute(:,:) ! is (nproma, blocks )
    TYPE(t_patch),POINTER            :: patch_2d
    CHARACTER(filename_max) :: OutputFileName   !< file name for reading in

    TYPE(t_stream_id) :: stream_id

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_testbed_read:ocean_test_read"

    patch_2d => patch_3d%p_patch_2d(1)

    !---------------------------------------------------------------------
    ! 2D Tests
    CALL openinputfile(stream_id, dynamics_grid_filename(1), patch_2d, &
      &                read_netcdf_broadcast_method)

    CALL read_2D(stream_id=stream_id, location=on_cells, &
      &          variable_name="cell_area", return_pointer=cell_data_broadcast)
    CALL read_2D(stream_id=stream_id, location=on_vertices, &
      &          variable_name="dual_area", return_pointer=vertex_data_broadcast)
    CALL read_2D(stream_id=stream_id, location=on_edges, &
      &          variable_name="phys_edge_id", return_pointer=edge_data_broadcast)

    CALL closeFile(stream_id)

    CALL openinputfile(stream_id, dynamics_grid_filename(1), patch_2d, &
      &                read_netcdf_distribute_method)

    CALL read_2D(stream_id=stream_id, location=on_cells, &
      &          variable_name="cell_area", return_pointer=cell_data_distribute)
    CALL read_2D(stream_id=stream_id, location=on_vertices, &
      &          variable_name="dual_area", return_pointer=vertex_data_distribute)
    CALL read_2D(stream_id=stream_id, location=on_edges, &
      &          variable_name="phys_edge_id", return_pointer=edge_data_distribute)

    CALL closeFile(stream_id)

    IF ( MAXVAL(ABS(cell_data_distribute - cell_data_broadcast )) > 0.0_wp ) &
      CALL finish(method_name, "Cell data check failed")
    IF ( MAXVAL(ABS(vertex_data_distribute - vertex_data_broadcast )) > 0.0_wp ) &
      CALL finish(method_name, "Vertex data check failed")
    IF ( MAXVAL(ABS(edge_data_distribute - edge_data_broadcast )) > 0.0_wp ) &
      CALL finish(method_name, "Edge data check failed")

    !---------------------------------------------------------------------
    CALL read_interface_test(read_netcdf_distribute_method)
    CALL read_interface_test(read_netcdf_broadcast_method)
    !---------------------------------------------------------------------

    RETURN ! disable the rest of tests (write/read) in order to run on buildbot
    !---------------------------------------------------------------------

    CALL openinputfile(stream_id, initialState_InputFileName, patch_2d, &
      &                read_netcdf_broadcast_method)

    CALL read_3D_time( stream_id=stream_id, location=on_cells, &
      &                variable_name="T", return_pointer=T )

    CALL closeFile(stream_id)
    !---------------------------------------------------------------------
    ! check

!     OutputFileName="testOut.nc"
!     CALL message(method_name,   OutputFileName)
!     return_status = netcdf_write_oncells_3d_time(  &
!       & filename=OutputFileName,                  &
!       & variable_name="T",                         &
!       & write_array=T,                             &
!       & patch=patch_2d )
!
!     CALL openinputfile(stream_id, OutputFileName, patch_2d, &
!       &                read_netcdf_broadcast_method)
!
!     CALL read_3D_time( stream_id=stream_id, location=on_cells, &
!       &                variable_name="T", return_pointer=T_check )
!     IF ( MAXVAL(ABS(T - T_check )) > 0.0_wp ) &
!       CALL finish(method_name, "Check failed")
!
!     CALL closeFile(stream_id)
    !---------------------------------------------------------------------


    DEALLOCATE(T)
    DEALLOCATE(T_check)

  CONTAINS

    SUBROUTINE read_interface_test(read_method)

      INTEGER, INTENT(IN) :: read_method

      TYPE(t_stream_id) :: stream_id
      TYPE(t_patch), TARGET :: dummy_patch
      INTEGER :: i
      REAL(wp), POINTER :: real_2d(:,:), real_3d(:,:,:), real_4d(:,:,:,:)
      INTEGER, POINTER :: int_2d(:,:), int_3d(:,:,:)

      CHARACTER(*), PARAMETER :: method_name = &
        "mo_ocean_testbed_read:read_interface_test"

      ! setup dummy patch
      dummy_patch%n_patch_cells_g = p_n_work * 2*nproma
      dummy_patch%n_patch_edges_g = p_n_work * 2*nproma
      dummy_patch%n_patch_verts_g = p_n_work * 2*nproma

      ALLOCATE(dummy_patch%cells%decomp_info%glb_index(2*nproma), &
        &      dummy_patch%edges%decomp_info%glb_index(2*nproma), &
        &      dummy_patch%verts%decomp_info%glb_index(2*nproma))

      dummy_patch%cells%decomp_info%glb_index(:) = &
        (/(i, i=p_pe_work+1, p_n_work * 2 * nproma, p_n_work)/)
      dummy_patch%edges%decomp_info%glb_index(:) = &
        (/(i, i=p_pe_work+1, p_n_work * 2 * nproma, p_n_work)/)
      dummy_patch%verts%decomp_info%glb_index(:) = &
        (/(i, i=p_pe_work+1, p_n_work * 2 * nproma, p_n_work)/)

      ALLOCATE(dummy_patch%cells%dist_io_data, dummy_patch%edges%dist_io_data, &
        &      dummy_patch%verts%dist_io_data)

      CALL setup_distrib_read(p_n_work * 2 * nproma, &
        &                     dummy_patch%cells%decomp_info, &
        &                     dummy_patch%cells%dist_io_data)
      CALL setup_distrib_read(p_n_work * 2 * nproma, &
        &                     dummy_patch%edges%decomp_info, &
        &                     dummy_patch%edges%dist_io_data)
      CALL setup_distrib_read(p_n_work * 2 * nproma, &
        &                     dummy_patch%verts%decomp_info, &
        &                     dummy_patch%verts%dist_io_data)

      dummy_patch%n_patch_cells = 2*nproma
      dummy_patch%n_patch_edges = 2*nproma
      dummy_patch%n_patch_verts = 2*nproma

      dummy_patch%comm_pat_scatter_c => &
        makeScatterPattern(1, dummy_patch%n_patch_cells, &
        &                  dummy_patch%cells%decomp_info%glb_index, &
        &                  p_comm_work)
      dummy_patch%comm_pat_scatter_e => &
        makeScatterPattern(1, dummy_patch%n_patch_edges, &
        &                  dummy_patch%edges%decomp_info%glb_index, &
        &                  p_comm_work)
      dummy_patch%comm_pat_scatter_v => &
        makeScatterPattern(1, dummy_patch%n_patch_verts, &
        &                  dummy_patch%verts%decomp_info%glb_index, &
        &                  p_comm_work)

      IF (p_pe_work == 0) &
        CALL generate_test_file("testfile.nc", nlev=10, ntime=5)

      CALL work_mpi_barrier()

      ! open input file
      CALL openinputfile(stream_id, "testfile.nc", dummy_patch, &
        &                read_method)

      CALL read_2D(stream_id, on_cells, 'cell_2d_real', return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2 * nproma) &
        CALL finish(method_name, "cell 2d real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2 * nproma, p_n_work)/))) &
        CALL finish(method_name, "cell 2d real return_point test failed")
      real_2d = 0
      CALL read_2D(stream_id, on_cells, 'cell_2d_real', fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "cell 2d real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test vertex 2d real
      CALL read_2D(stream_id, on_vertices, 'vertex_2d_real', return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, "vertex 2d real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d real return_point test failed")
      real_2d = 0
      CALL read_2D(stream_id, on_vertices, 'vertex_2d_real', fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test edge 2d real
      CALL read_2D(stream_id, on_edges, 'edge_2d_real', return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, "edge 2d real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "edge 2d real return_point test failed")
      real_2d = 0
      CALL read_2D(stream_id, on_edges, 'edge_2d_real', fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "edge 2d real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test cell 2d int
      CALL read_2D_int(stream_id, on_cells, 'cell_2d_int', return_pointer=int_2d)
      IF (SIZE(int_2d) /= 2*nproma) &
        CALL finish(method_name, "cell 2d int return_pointer size test failed")
      IF (ANY((/(int_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "cell 2d int return_point test failed")
      int_2d = 0
      CALL read_2D_int(stream_id, on_cells, 'cell_2d_int', fill_array=int_2d)
      IF (ANY((/(int_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "cell 2d int fill_array test failed")
      DEALLOCATE(int_2d)

      ! test vertex 2d int
      CALL read_2D_int(stream_id, on_vertices, 'vertex_2d_int', return_pointer=int_2d)
      IF (SIZE(int_2d) /= 2*nproma) &
        CALL finish(method_name, "vertex 2d int return_pointer size test failed")
      IF (ANY((/(int_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d int return_point test failed")
      int_2d = 0
      CALL read_2D_int(stream_id, on_vertices, 'vertex_2d_int', fill_array=int_2d)
      IF (ANY((/(int_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d int fill_array test failed")
      DEALLOCATE(int_2d)

      ! test edge 2d int
      CALL read_2D_int(stream_id, on_edges, 'edge_2d_int', return_pointer=int_2d)
      IF (SIZE(int_2d) /= 2*nproma) &
        CALL finish(method_name, "edge 2d int return_pointer size test failed")
      IF (ANY((/(int_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "edge 2d int return_point test failed")
      int_2d = 0
      CALL read_2D_int(stream_id, on_edges, 'edge_2d_int', fill_array=int_2d)
      IF (ANY((/(int_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "edge 2d int fill_array test failed")
      DEALLOCATE(int_2d)

      ! test cell 3d real
      CALL read_3d(stream_id, on_cells, 'cell_3d_real', return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 10) &
        CALL finish(method_name, "cell 3d real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 10, p_n_work)/))) &
        CALL finish(method_name, "cell 3d real return_point test failed")
      real_3d = 0
      CALL read_3d(stream_id, on_cells, 'cell_3d_real', fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 10, p_n_work)/))) &
        CALL finish(method_name, "cell 3d real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test vertex 3d real
      CALL read_3d(stream_id, on_vertices, 'vertex_3d_real', return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 10) &
        CALL finish(method_name, "vertex 3d real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 10, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d real return_point test failed")
      real_3d = 0
      CALL read_3d(stream_id, on_vertices, 'vertex_3d_real', fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 10, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test edge 3d real
      CALL read_3d(stream_id, on_edges, 'edge_3d_real', return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 10) &
        CALL finish(method_name, "edge 3d real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 10, p_n_work)/))) &
        CALL finish(method_name, "edge 3d real return_point test failed")
      real_3d = 0
      CALL read_3d(stream_id, on_edges, 'edge_3d_real', fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 10, p_n_work)/))) &
        CALL finish(method_name, "edge 3d real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test cell 2d time real
      CALL read_2d_time(stream_id, on_cells, 'cell_2d_time_real', &
        &               return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "cell 2d time real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 2d time real return_point test failed")
      real_3d = 0
      CALL read_2d_time(stream_id, on_cells, 'cell_2d_time_real', &
        &               fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                           (i-1)/(p_n_work*2*nproma)+1) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 2d time real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test vertex 2d time real
      CALL read_2d_time(stream_id, on_vertices, 'vertex_2d_time_real', &
        &               return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "vertex 2d time real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d time real return_point test failed")
      real_3d = 0
      CALL read_2d_time(stream_id, on_vertices, 'vertex_2d_time_real', &
        &               fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d time real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test edge 2d time real
      CALL read_2d_time(stream_id, on_edges, 'edge_2d_time_real', &
        &               return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "edge 2d time real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 2d time real return_point test failed")
      real_3d = 0
      CALL read_2d_time(stream_id, on_edges, 'edge_2d_time_real', &
        &               fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 2d time real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test cell 2d time real (section)
      CALL read_2d_time(stream_id, on_cells, 'cell_2d_time_real', &
        &               return_pointer=real_3d, start_timestep=2, end_timestep=4)
      IF (SIZE(real_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "cell 2d time real (section) return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 2d time real (section) return_point test failed")
      real_3d = 0
      CALL read_2d_time(stream_id, on_cells, 'cell_2d_time_real', &
        &               fill_array=real_3d, start_timestep=2, end_timestep=4)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 2d time real (section) fill_array test failed")
      DEALLOCATE(real_3d)

      ! test vertex 2d time real (section)
      CALL read_2d_time(stream_id, on_vertices, 'vertex_2d_time_real', &
        &               return_pointer=real_3d, start_timestep=2, end_timestep=4)
      IF (SIZE(real_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "vertex 2d time real (section) return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 2d time real (section) return_point test failed")
      real_3d = 0
      CALL read_2d_time(stream_id, on_vertices, 'vertex_2d_time_real', &
        &               fill_array=real_3d, start_timestep=2, end_timestep=4)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 2d time real (section) fill_array test failed")
      DEALLOCATE(real_3d)

      ! test edge 2d time real (section)
      CALL read_2d_time(stream_id, on_edges, 'edge_2d_time_real', &
        &               return_pointer=real_3d, start_timestep=2, end_timestep=4)
      IF (SIZE(real_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "edge 2d time real (section) return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 2d time real (section) return_point test failed")
      real_3d = 0
      CALL read_2d_time(stream_id, on_edges, 'edge_2d_time_real', &
        &               fill_array=real_3d, start_timestep=2, end_timestep=4)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 2d time real (section) fill_array test failed")
      DEALLOCATE(real_3d)

      ! test cell 2d extdim real
      CALL read_2d_extdim(stream_id, on_cells, 'cell_2d_time_real', &
        &                 return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "cell 2d extdim real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 2d extdim real return_point test failed")
      real_3d = 0
      CALL read_2d_extdim(stream_id, on_cells, 'cell_2d_time_real', &
        &                 fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 2d extdim real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test vertex 2d extdim real
      CALL read_2d_extdim(stream_id, on_vertices, 'vertex_2d_time_real', &
        &                 return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "vertex 2d extdim real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d extdim real return_point test failed")
      real_3d = 0
      CALL read_2d_extdim(stream_id, on_vertices, 'vertex_2d_time_real', &
        &                 fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d extdim real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test edge 2d extdim real
      CALL read_2d_extdim(stream_id, on_edges, 'edge_2d_time_real', &
        &                 return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "edge 2d extdim real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 2d extdim real return_point test failed")
      real_3d = 0
      CALL read_2d_extdim(stream_id, on_edges, 'edge_2d_time_real', &
        &                 fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 2d extdim real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test cell 2d extdim real (section)
      CALL read_2d_extdim(stream_id, on_cells, 'cell_2d_time_real', &
        &                 return_pointer=real_3d, start_extdim=2, &
        &                 end_extdim=4, extdim_name='time')
      IF (SIZE(real_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "cell 2d extdim real (section) return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 2d time real (section) return_point test failed")
      real_3d = 0
      CALL read_2d_extdim(stream_id, on_cells, 'cell_2d_time_real', &
        &                 fill_array=real_3d, start_extdim=2, end_extdim=4, &
        &                 extdim_name='time')
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 2d extdim real (section) fill_array test failed")
      DEALLOCATE(real_3d)

      ! test vertex 2d extdim real (section)
      CALL read_2d_extdim(stream_id, on_vertices, 'vertex_2d_time_real', &
        &                 return_pointer=real_3d, start_extdim=2, &
        &                 end_extdim=4, extdim_name='time')
      IF (SIZE(real_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "vertex 2d extdim real (section) return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 2d time real (section) return_point test failed")
      real_3d = 0
      CALL read_2d_extdim(stream_id, on_vertices, 'vertex_2d_time_real', &
        &                 fill_array=real_3d, start_extdim=2, end_extdim=4, &
        &                 extdim_name='time')
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 2d extdim real (section) fill_array test failed")
      DEALLOCATE(real_3d)

      ! test edge 2d extdim real (section)
      CALL read_2d_extdim(stream_id, on_edges, 'edge_2d_time_real', &
        &                 return_pointer=real_3d, start_extdim=2, &
        &                 end_extdim=4, extdim_name='time')
      IF (SIZE(real_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "edge 2d extdim real (section) return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 2d time real (section) return_point test failed")
      real_3d = 0
      CALL read_2d_extdim(stream_id, on_edges, 'edge_2d_time_real', &
        &                 fill_array=real_3d, start_extdim=2, end_extdim=4, &
        &                 extdim_name='time')
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 2d extdim real (section) fill_array test failed")
      DEALLOCATE(real_3d)

      ! test cell 2d extdim int
      CALL read_2d_extdim_int(stream_id, on_cells, 'cell_2d_time_int', &
        &                     return_pointer=int_3d)
      IF (SIZE(int_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "cell 2d extdim int return_pointer size test failed")
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 2d extdim int return_point test failed")
      int_3d = 0
      CALL read_2d_extdim_int(stream_id, on_cells, 'cell_2d_time_int', &
        &                     fill_array=int_3d)
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 2d extdim int fill_array test failed")
      DEALLOCATE(int_3d)

      ! test vertex 2d extdim int
      CALL read_2d_extdim_int(stream_id, on_vertices, 'vertex_2d_time_int', &
        &                     return_pointer=int_3d)
      IF (SIZE(int_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "vertex 2d extdim int return_pointer size test failed")
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &        blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &        (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d extdim int return_point test failed")
      int_3d = 0
      CALL read_2d_extdim_int(stream_id, on_vertices, 'vertex_2d_time_int', &
        &                     fill_array=int_3d)
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d extdim int fill_array test failed")
      DEALLOCATE(int_3d)

      ! test edge 2d extdim int
      CALL read_2d_extdim_int(stream_id, on_edges, 'edge_2d_time_int', &
        &                     return_pointer=int_3d)
      IF (SIZE(int_3d) /= 2*nproma * 5) &
        CALL finish(method_name, &
          &         "edge 2d extdim int return_pointer size test failed")
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 2d extdim int return_point test failed")
      int_3d = 0
      CALL read_2d_extdim_int(stream_id, on_edges, 'edge_2d_time_int', &
        &                     fill_array=int_3d)
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 2d extdim int fill_array test failed")
      DEALLOCATE(int_3d)

      ! test cell 2d extdim int (section)
      CALL read_2d_extdim_int(stream_id, on_cells, 'cell_2d_time_int', &
        &                     return_pointer=int_3d, start_extdim=2, &
        &                 end_extdim=4, extdim_name='time')
      IF (SIZE(int_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "cell 2d extdim int (section) return_pointer size test failed")
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 2d time int (section) return_point test failed")
      int_3d = 0
      CALL read_2d_extdim_int(stream_id, on_cells, 'cell_2d_time_int', &
        &                     fill_array=int_3d, start_extdim=2, end_extdim=4, &
        &                     extdim_name='time')
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 2d extdim int (section) fill_array test failed")
      DEALLOCATE(int_3d)

      ! test vertex 2d extdim int (section)
      CALL read_2d_extdim_int(stream_id, on_vertices, 'vertex_2d_time_int', &
        &                     return_pointer=int_3d, start_extdim=2, &
        &                     end_extdim=4, extdim_name='time')
      IF (SIZE(int_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "vertex 2d extdim int (section) return_pointer size test failed")
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 2d time int (section) return_point test failed")
      int_3d = 0
      CALL read_2d_extdim_int(stream_id, on_vertices, 'vertex_2d_time_int', &
        &                     fill_array=int_3d, start_extdim=2, end_extdim=4, &
        &                     extdim_name='time')
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 2d extdim int (section) fill_array test failed")
      DEALLOCATE(int_3d)

      ! test edge 2d extdim int (section)
      CALL read_2d_extdim_int(stream_id, on_edges, 'edge_2d_time_int', &
        &                     return_pointer=int_3d, start_extdim=2, &
        &                     end_extdim=4, extdim_name='time')
      IF (SIZE(int_3d) /= 2*nproma * 3) &
        CALL finish(method_name, &
          &         "edge 2d extdim int (section) return_pointer size test failed")
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 2d time int (section) return_point test failed")
      int_3d = 0
      CALL read_2d_extdim_int(stream_id, on_edges, 'edge_2d_time_int', &
        &                     fill_array=int_3d, start_extdim=2, end_extdim=4, &
        &                     extdim_name='time')
      IF (ANY((/(int_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma)+1) /= i + p_n_work*2*nproma, &
        &        i=p_pe_work+1, p_n_work * 2*nproma * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 2d extdim int (section) fill_array test failed")
      DEALLOCATE(int_3d)

      ! test cell 3d time real
      CALL read_3d_time(stream_id, on_cells, 'cell_3d_time_real', &
        &               return_pointer=real_4d)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 5) &
        CALL finish(method_name, &
          &         "cell 3d time real return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 3d time real return_point test failed")
      real_4d = 0
      CALL read_3d_time(stream_id, on_cells, 'cell_3d_time_real', &
        &               fill_array=real_4d)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 3d time real fill_array test failed")
      DEALLOCATE(real_4d)

      ! test vertex 3d time real
      CALL read_3d_time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &               return_pointer=real_4d)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 5) &
        CALL finish(method_name, &
          &         "vertex 3d time real return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d time real return_point test failed")
      real_4d = 0
      CALL read_3d_time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &               fill_array=real_4d)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d time real fill_array test failed")
      DEALLOCATE(real_4d)

      ! test edge 3d time real
      CALL read_3d_time(stream_id, on_edges, 'edge_3d_time_real', &
        &               return_pointer=real_4d)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 5) &
        CALL finish(method_name, &
          &         "edge 3d time real return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 3d time real return_point test failed")
      real_4d = 0
      CALL read_3d_time(stream_id, on_edges, 'edge_3d_time_real', &
        &               fill_array=real_4d)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 3d time real fill_array test failed")
      DEALLOCATE(real_4d)

      ! test cell 3d time real (section)
      CALL read_3d_time(stream_id, on_cells, 'cell_3d_time_real', &
        &               return_pointer=real_4d, start_timestep=2, end_timestep=4)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 3) &
        CALL finish(method_name, &
          &         "cell 3d time real (section) return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) THEN
        CALL finish(method_name, &
          &         "cell 3d time real (section) return_point test failed")
      END IF
      real_4d = 0
      CALL read_3d_time(stream_id, on_cells, 'cell_3d_time_real', &
        &               fill_array=real_4d, start_timestep=2, end_timestep=4)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 3d time real (section) fill_array test failed")
      DEALLOCATE(real_4d)

      ! test vertex 3d time real (section)
      CALL read_3d_time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &               return_pointer=real_4d, start_timestep=2, end_timestep=4)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 3) &
        CALL finish(method_name, &
          &         "vertex 3d time real (section) return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) THEN
        CALL finish(method_name, &
          &         "vertex 3d time real (section) return_point test failed")
      END IF
      real_4d = 0
      CALL read_3d_time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &               fill_array=real_4d, start_timestep=2, end_timestep=4)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 3d time real (section) fill_array test failed")
      DEALLOCATE(real_4d)

      ! test edge 3d time real (section)
      CALL read_3d_time(stream_id, on_edges, 'edge_3d_time_real', &
        &               return_pointer=real_4d, start_timestep=2, end_timestep=4)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 3) &
        CALL finish(method_name, &
          &         "edge 3d time real (section) return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) THEN
        CALL finish(method_name, &
          &         "edge 3d time real (section) return_point test failed")
      END IF
      real_4d = 0
      CALL read_3d_time(stream_id, on_edges, 'edge_3d_time_real', &
        &               fill_array=real_4d, start_timestep=2, end_timestep=4)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 3d time real (section) fill_array test failed")
      DEALLOCATE(real_4d)

      ! test cell 3d extdim real
      CALL read_3d_extdim(stream_id, on_cells, 'cell_3d_time_real', &
        &                 return_pointer=real_4d)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 5) &
        CALL finish(method_name, &
          &         "cell 3d extdim real return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
                &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 3d extdim real return_point test failed")
      real_4d = 0
      CALL read_3d_extdim(stream_id, on_cells, 'cell_3d_time_real', &
        &                 fill_array=real_4d)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
                &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "cell 3d extdim real fill_array test failed")
      DEALLOCATE(real_4d)

      ! test vertex 3d extdim real
      CALL read_3d_extdim(stream_id, on_vertices, 'vertex_3d_time_real', &
        &                 return_pointer=real_4d)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 5) &
        CALL finish(method_name, &
          &         "vertex 3d extdim real return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d extdim real return_point test failed")
      real_4d = 0
      CALL read_3d_extdim(stream_id, on_vertices, 'vertex_3d_time_real', &
        &                 fill_array=real_4d)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
                &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d extdim real fill_array test failed")
      DEALLOCATE(real_4d)

      ! test edge 3d extdim real
      CALL read_3d_extdim(stream_id, on_edges, 'edge_3d_time_real', &
        &                 return_pointer=real_4d)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 5) &
        CALL finish(method_name, &
          &         "edge 3d extdim real return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
                &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 3d extdim real return_point test failed")
      real_4d = 0
      CALL read_3d_extdim(stream_id, on_edges, 'edge_3d_time_real', &
        &                 fill_array=real_4d)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
                &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 5, p_n_work)/))) &
        CALL finish(method_name, "edge 3d extdim real fill_array test failed")
      DEALLOCATE(real_4d)

      ! test cell 3d extdim real (section)
      CALL read_3d_extdim(stream_id, on_cells, 'cell_3d_time_real', &
        &               return_pointer=real_4d, start_extdim=2, end_extdim=4)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 3) &
        CALL finish(method_name, &
          &         "cell 3d extdim real (section) return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) THEN
        CALL finish(method_name, &
          &         "cell 3d extdim real (section) return_point test failed")
      END IF
      real_4d = 0
      CALL read_3d_extdim(stream_id, on_cells, 'cell_3d_time_real', &
        &               fill_array=real_4d, start_extdim=2, end_extdim=4)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 3d extdim real (section) fill_array test failed")
      DEALLOCATE(real_4d)

      ! test vertex 3d extdim real (section)
      CALL read_3d_extdim(stream_id, on_vertices, 'vertex_3d_time_real', &
        &               return_pointer=real_4d, start_extdim=2, end_extdim=4)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 3) &
        CALL finish(method_name, &
          &         "vertex 3d extdim real (section) return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) THEN
        CALL finish(method_name, &
          &         "vertex 3d extdim real (section) return_point test failed")
      END IF
      real_4d = 0
      CALL read_3d_extdim(stream_id, on_vertices, 'vertex_3d_time_real', &
        &               fill_array=real_4d, start_extdim=2, end_extdim=4)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 3d extdim real (section) fill_array test failed")
      DEALLOCATE(real_4d)

      ! test edge 3d extdim real (section)
      CALL read_3d_extdim(stream_id, on_edges, 'edge_3d_time_real', &
        &               return_pointer=real_4d, start_extdim=2, end_extdim=4)
      IF (SIZE(real_4d) /= 2*nproma * 10 * 3) &
        CALL finish(method_name, &
          &         "edge 3d extdim real (section) return_pointer size test failed")
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &               (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) THEN
        CALL finish(method_name, &
          &         "edge 3d extdim real (section) return_point test failed")
      END IF
      real_4d = 0
      CALL read_3d_extdim(stream_id, on_edges, 'edge_3d_time_real', &
        &               fill_array=real_4d, start_extdim=2, end_extdim=4)
      IF (ANY((/(real_4d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                MOD((i-1)/(p_n_work*2*nproma),10)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                (i-1)/(p_n_work*2*nproma*10)+1) /= i + p_n_work*2*nproma*10, &
        &                i=p_pe_work+1, p_n_work * 2*nproma * 10 * 3, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 3d extdim real (section) fill_array test failed")
      DEALLOCATE(real_4d)

      ! test cell 3d extdim real with dimension names
      CALL read_3d_extdim(stream_id, on_cells, 'cell_3d_time_real', &
        &                 return_pointer=real_4d, levelsDimName='levels', &
        &                 extdim_name='time')
      IF (SIZE(real_4d) /= 2*nproma * 10 * 5) &
        CALL finish(method_name, &
          &         "cell 3d extdim real with dimension names return_pointer size test failed")
      DEALLOCATE(real_4d)

      ! close input file
      CALL closeFile(stream_id)

      CALL work_mpi_barrier()

      IF (p_pe_work == 0) THEN
        return_status_c = util_unlink("testfile.nc")
        CALL generate_test_file("testfile.nc", nlev=1, ntime=1)
      END IF

      CALL work_mpi_barrier()

      ! open input file
      CALL openinputfile(stream_id, "testfile.nc", dummy_patch, read_method)

      ! test cell 2d 1time real
      CALL read_2D_1time(stream_id, on_cells, 'cell_2d_time_real', return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, "cell 2d 1time real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "cell 2d real 1time return_point test failed")
      real_2d = 0
      CALL read_2D_1time(stream_id, on_cells, 'cell_2d_time_real', fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "cell 2d 1time real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test vertex 2d 1time real
      CALL read_2D_1time(stream_id, on_vertices, 'vertex_2d_time_real', &
        &                return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, "vertex 2d 1time real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d 1time real return_point test failed")
      real_2d = 0
      CALL read_2D_1time(stream_id, on_vertices, 'vertex_2d_time_real', &
        &                fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d 1time real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test edge 2d 1time real
      CALL read_2D_1time(stream_id, on_edges, 'edge_2d_time_real', &
        &                return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, "edge 2d 1time real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "edge 2d 1time real return_point test failed")
      real_2d = 0
      CALL read_2D_1time(stream_id, on_edges, 'edge_2d_time_real', fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "edge 2d 1time real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test cell 2d 1lev 1time real
      CALL read_2D_1lev_1time(stream_id, on_cells, 'cell_3d_time_real', &
        &                     return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, &
          &         "cell 2d 1lev 1time real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, &
          &         "cell 2d 1lev 1time real return_point test failed")
      real_2d = 0
      CALL read_2D_1lev_1time(stream_id, on_cells, 'cell_3d_time_real', &
        &                     fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "cell 2d 1lev 1time real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test vertex 2d 1lev 1time real
      CALL read_2D_1lev_1time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &                     return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, &
          &         "vertex 2d 1lev 1time real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, &
          &         "vertex 2d 1lev 1time real return_point test failed")
      real_2d = 0
      CALL read_2D_1lev_1time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &                     fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "vertex 2d 1lev 1time real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test edge 2d 1lev 1time real
      CALL read_2D_1lev_1time(stream_id, on_edges, 'edge_3d_time_real', &
        &                     return_pointer=real_2d)
      IF (SIZE(real_2d) /= 2*nproma) &
        CALL finish(method_name, &
          &         "edge 2d 1lev 1time real return_pointer size test failed")
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, &
          &         "edge 2d 1lev 1time real return_point test failed")
      real_2d = 0
      CALL read_2D_1lev_1time(stream_id, on_edges, 'edge_3d_time_real', &
        &                     fill_array=real_2d)
      IF (ANY((/(real_2d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, &
        &        i=p_pe_work+1, p_n_work * 2*nproma, p_n_work)/))) &
        CALL finish(method_name, "edge 2d 1lev 1time real fill_array test failed")
      DEALLOCATE(real_2d)

      ! test cell 3d 1time real
      CALL read_3d_1time(stream_id, on_cells, 'cell_3d_time_real', &
        &                return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 1) &
        CALL finish(method_name, &
          &         "cell 3d 1time real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 1, p_n_work)/))) &
        CALL finish(method_name, "cell 3d 1time real return_point test failed")
      real_3d = 0
      CALL read_3d_1time(stream_id, on_cells, 'cell_3d_time_real', &
        &                fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 1, p_n_work)/))) &
        CALL finish(method_name, "cell 3d 1time real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test vertex 3d 1time real
      CALL read_3d_1time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &                return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 1) &
        CALL finish(method_name, &
          &         "vertex 3d 1time real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 1, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d 1time real return_point test failed")
      real_3d = 0
      CALL read_3d_1time(stream_id, on_vertices, 'vertex_3d_time_real', &
        &                fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 1, p_n_work)/))) &
        CALL finish(method_name, "vertex 3d 1time real fill_array test failed")
      DEALLOCATE(real_3d)

      ! test edge 3d 1time real
      CALL read_3d_1time(stream_id, on_edges, 'edge_3d_time_real', &
        &                return_pointer=real_3d)
      IF (SIZE(real_3d) /= 2*nproma * 1) &
        CALL finish(method_name, &
          &         "edge 3d 1time real return_pointer size test failed")
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 1, p_n_work)/))) &
        CALL finish(method_name, "edge 3d 1time real return_point test failed")
      real_3d = 0
      CALL read_3d_1time(stream_id, on_edges, 'edge_3d_time_real', &
        &                fill_array=real_3d)
      IF (ANY((/(real_3d(idx_no(MOD((i-1)/p_n_work,2*nproma)+1), (i-1)/(p_n_work*2*nproma)+1, &
        &                blk_no(MOD((i-1)/p_n_work,2*nproma)+1)) /= i, i=p_pe_work+1, &
        &        p_n_work * 2*nproma * 1, p_n_work)/))) &
        CALL finish(method_name, "edge 3d 1time real fill_array test failed")
      DEALLOCATE(real_3d)

      ! close input file
      CALL closeFile(stream_id)

      CALL work_mpi_barrier()

      IF (p_pe_work == 0)  return_status_c = util_unlink("testfile.nc")

      ! clean up dummy patch
      CALL delete_distrib_read(dummy_patch%cells%dist_io_data)
      CALL delete_distrib_read(dummy_patch%edges%dist_io_data)
      CALL delete_distrib_read(dummy_patch%verts%dist_io_data)

      DEALLOCATE(dummy_patch%cells%dist_io_data, &
        &        dummy_patch%edges%dist_io_data, &
        &        dummy_patch%verts%dist_io_data)

      DEALLOCATE(dummy_patch%cells%decomp_info%glb_index, &
        &        dummy_patch%edges%decomp_info%glb_index, &
        &        dummy_patch%verts%decomp_info%glb_index)

      CALL dummy_patch%comm_pat_scatter_c%destruct()
      CALL dummy_patch%comm_pat_scatter_e%destruct()
      CALL dummy_patch%comm_pat_scatter_v%destruct()
    END SUBROUTINE read_interface_test

    SUBROUTINE generate_test_file(filename, nlev, ntime)

      CHARACTER(len=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: nlev, ntime

      CHARACTER(*), PARAMETER :: method_name = &
        "mo_ocean_testbed_read:generate_test_file"

      INTEGER :: ncid, i

      INTEGER :: cell_dimid, vertex_dimid, edge_dimid, level_dimid, time_dimid

      ! We will write surface temperature and pressure fields.
      INTEGER :: cell_2d_real_varid, vertex_2d_real_varid, &
       &         edge_2d_real_varid
      INTEGER :: cell_3d_real_varid, vertex_3d_real_varid, &
       &         edge_3d_real_varid
      INTEGER :: cell_2d_time_real_varid, vertex_2d_time_real_varid, &
       &         edge_2d_time_real_varid
      INTEGER :: cell_3d_time_real_varid, vertex_3d_time_real_varid, &
       &         edge_3d_time_real_varid
      INTEGER :: cell_2d_int_varid, vertex_2d_int_varid, &
       &         edge_2d_int_varid
      INTEGER :: cell_3d_int_varid, vertex_3d_int_varid, &
       &         edge_3d_int_varid
      INTEGER :: cell_2d_time_int_varid, vertex_2d_time_int_varid, &
       &         edge_2d_time_int_varid
      INTEGER :: cell_3d_time_int_varid, vertex_3d_time_int_varid, &
       &         edge_3d_time_int_varid

      ! create the file.
      call nf( nf_create(filename, NF_CLOBBER, ncid), method_name )

      ! define the dimensions.
      CALL nf(nf_def_dim(ncid, 'cell', p_n_work * 2*nproma, cell_dimid), &
       &      method_name)
      CALL nf(nf_def_dim(ncid, 'vertex', p_n_work * 2*nproma, vertex_dimid), &
       &      method_name)
      CALL nf(nf_def_dim(ncid, 'edge', p_n_work * 2*nproma, edge_dimid), &
       &      method_name)
      CALL nf(nf_def_dim(ncid, 'levels', nlev, level_dimid), &
       &      method_name)
      CALL nf(nf_def_dim(ncid, 'time', ntime, time_dimid), &
       &      method_name)

      ! define variables

      CALL nf(nf_def_var(ncid, 'cell_2d_real', NF_DOUBLE, 1, (/cell_dimid/), &
        &                cell_2d_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_2d_real', NF_DOUBLE, 1, (/vertex_dimid/),&
        &                vertex_2d_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_2d_real', NF_DOUBLE, 1, (/edge_dimid/), &
        &                edge_2d_real_varid), method_name )

      CALL nf(nf_def_var(ncid, 'cell_3d_real', NF_DOUBLE, 2, &
        &                (/cell_dimid, level_dimid/), &
        &                cell_3d_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_3d_real', NF_DOUBLE, 2, &
        &                (/vertex_dimid, level_dimid/), &
        &                vertex_3d_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_3d_real', NF_DOUBLE, 2, &
        &                (/edge_dimid, level_dimid/), &
        &                edge_3d_real_varid), method_name )

      CALL nf(nf_def_var(ncid, 'cell_2d_time_real', NF_DOUBLE, 2, &
        &                (/cell_dimid, time_dimid/), &
        &                cell_2d_time_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_2d_time_real', NF_DOUBLE, 2, &
        &                (/vertex_dimid, time_dimid/), &
        &                vertex_2d_time_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_2d_time_real', NF_DOUBLE, 2, &
        &                (/edge_dimid, time_dimid/), &
        &                edge_2d_time_real_varid), method_name )

      CALL nf(nf_def_var(ncid, 'cell_3d_time_real', NF_DOUBLE, 3, &
        &                (/cell_dimid, level_dimid, time_dimid/), &
        &                cell_3d_time_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_3d_time_real', NF_DOUBLE, 3, &
        &                (/vertex_dimid, level_dimid, time_dimid/), &
        &                vertex_3d_time_real_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_3d_time_real', NF_DOUBLE, 3, &
        &                (/edge_dimid, level_dimid, time_dimid/), &
        &                edge_3d_time_real_varid), method_name )

      CALL nf(nf_def_var(ncid, 'cell_2d_int', NF_INT, 1, (/cell_dimid/), &
        &                cell_2d_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_2d_int', NF_INT, 1, (/vertex_dimid/), &
        &                vertex_2d_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_2d_int', NF_INT, 1, (/edge_dimid/), &
        &                edge_2d_int_varid), method_name )

      CALL nf(nf_def_var(ncid, 'cell_3d_int', NF_INT, 2, &
        &                (/cell_dimid, level_dimid/), &
        &                cell_3d_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_3d_int', NF_INT, 2, &
        &                (/vertex_dimid, level_dimid/), &
        &                vertex_3d_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_3d_int', NF_INT, 2, &
        &                (/edge_dimid, level_dimid/), &
        &                edge_3d_int_varid), method_name )

      CALL nf(nf_def_var(ncid, 'cell_2d_time_int', NF_INT, 2, &
        &                (/cell_dimid, time_dimid/), &
        &                cell_2d_time_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_2d_time_int', NF_INT, 2, &
        &                (/vertex_dimid, time_dimid/), &
        &                vertex_2d_time_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_2d_time_int', NF_INT, 2, &
        &                (/edge_dimid, time_dimid/), &
        &                edge_2d_time_int_varid), method_name )

      CALL nf(nf_def_var(ncid, 'cell_3d_time_int', NF_INT, 3, &
        &                (/cell_dimid, level_dimid, time_dimid/), &
        &                cell_3d_time_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'vertex_3d_time_int', NF_INT, 3, &
        &                (/vertex_dimid, level_dimid, time_dimid/), &
        &                vertex_3d_time_int_varid), method_name )
      CALL nf(nf_def_var(ncid, 'edge_3d_time_int', NF_INT, 3, &
        &                (/edge_dimid, level_dimid, time_dimid/), &
        &                edge_3d_time_int_varid), method_name )

       ! end define mode
       call nf( nf_enddef(ncid), method_name )

       ! fill in dummy values

       CALL nf(nf_put_var_double(ncid, cell_2d_real_varid, &
        &                        (/(REAL(i,wp),i=1,p_n_work * 2*nproma)/)), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, vertex_2d_real_varid, &
        &                        (/(REAL(i,wp),i=1,p_n_work * 2*nproma)/)), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, edge_2d_real_varid, &
        &                        (/(REAL(i,wp),i=1,p_n_work * 2*nproma)/)), &
        &      method_name)

       CALL nf(nf_put_var_double(ncid, cell_3d_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * nlev)/), &
        &                                (/p_n_work * 2*nproma, nlev/))), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, vertex_3d_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * nlev)/), &
        &                                (/p_n_work * 2*nproma, nlev/))), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, edge_3d_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * nlev)/), &
        &                                (/p_n_work * 2*nproma, nlev/))), &
        &      method_name)

       CALL nf(nf_put_var_double(ncid, cell_2d_time_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * ntime)/), &
        &                                (/p_n_work * 2*nproma, ntime/))), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, vertex_2d_time_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * ntime)/), &
        &                                (/p_n_work * 2*nproma, ntime/))), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, edge_2d_time_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * ntime)/), &
        &                                (/p_n_work * 2*nproma, ntime/))), &
        &      method_name)

       CALL nf(nf_put_var_double(ncid, cell_3d_time_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * nlev * ntime)/), &
        &                                (/p_n_work * 2*nproma, nlev, ntime/))), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, vertex_3d_time_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * nlev * ntime)/), &
        &                                (/p_n_work * 2*nproma, nlev, ntime/))), &
        &      method_name)
       CALL nf(nf_put_var_double(ncid, edge_3d_time_real_varid, &
        &                        RESHAPE((/(REAL(i,wp),i=1,p_n_work * 2 * &
        &                                   nproma * nlev * ntime)/), &
        &                                (/p_n_work * 2*nproma, nlev, ntime/))), &
        &      method_name)

       CALL nf(nf_put_var_int(ncid, cell_2d_int_varid, &
        &                     (/(i,i=1,p_n_work * 2*nproma)/)), method_name)
       CALL nf(nf_put_var_int(ncid, vertex_2d_int_varid, &
        &                     (/(i,i=1,p_n_work * 2*nproma)/)), method_name)
       CALL nf(nf_put_var_int(ncid, edge_2d_int_varid, &
        &                     (/(i,i=1,p_n_work * 2*nproma)/)), method_name)

       CALL nf(nf_put_var_int(ncid, cell_3d_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2*nproma * nlev)/), &
        &                           (/p_n_work * 2*nproma, nlev/))), method_name)
       CALL nf(nf_put_var_int(ncid, vertex_3d_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2*nproma * nlev)/), &
        &                           (/p_n_work * 2*nproma, nlev/))), method_name)
       CALL nf(nf_put_var_int(ncid, edge_3d_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2*nproma * nlev)/), &
        &                           (/p_n_work * 2*nproma, nlev/))), method_name)

       CALL nf(nf_put_var_int(ncid, cell_2d_time_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2*nproma * ntime)/), &
        &                           (/p_n_work * 2*nproma, ntime/))), method_name)
       CALL nf(nf_put_var_int(ncid, vertex_2d_time_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2*nproma * ntime)/), &
        &                           (/p_n_work * 2*nproma, ntime/))), method_name)
       CALL nf(nf_put_var_int(ncid, edge_2d_time_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2*nproma * ntime)/), &
        &                           (/p_n_work * 2*nproma, ntime/))), method_name)

       CALL nf(nf_put_var_int(ncid, cell_3d_time_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2 * nproma * nlev * ntime)/),&
        &                             (/p_n_work * 2*nproma, nlev, ntime/))), &
        &      method_name)
       CALL nf(nf_put_var_int(ncid, vertex_3d_time_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2 * nproma * nlev * ntime)/),&
        &                             (/p_n_work * 2*nproma, nlev, ntime/))), &
        &      method_name)
       CALL nf(nf_put_var_int(ncid, edge_3d_time_int_varid, &
        &                     RESHAPE((/(i,i=1,p_n_work * 2 * nproma * nlev * ntime)/),&
        &                             (/p_n_work * 2*nproma, nlev, ntime/))), &
        &      method_name)

       ! close the file
       call nf( nf_close(ncid), method_name )

    END SUBROUTINE

  END SUBROUTINE ocean_test_read
  !-------------------------------------------------------------------------

END MODULE mo_ocean_testbed_read

