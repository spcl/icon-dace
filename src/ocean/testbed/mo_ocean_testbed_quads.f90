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

! @brief tests the mo_read_interface methods

MODULE mo_ocean_testbed_quads

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
  USE mo_read_interface
  USE mo_read_netcdf_distributed, ONLY: setup_distrib_read, delete_distrib_read
  USE mo_communication,       ONLY: idx_no, blk_no, makeScatterPattern
  USE mo_util_file,           ONLY: util_unlink

!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

  PUBLIC :: ocean_test_quads

CONTAINS

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE ocean_test_quads(namelist_filename, shr_namelist_filename, patch_3d)

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

    CHARACTER(*), PARAMETER :: method_name = "mo_ocean_testbed_quads:ocean_test_quads"

    patch_2d => patch_3d%p_patch_2d(1)

    !---------------------------------------------------------------------
    ! 2D Tests
    CALL openinputfile(stream_id, dynamics_grid_filename(1), patch_2d, &
      &                read_netcdf_broadcast_method)

    CALL read_2D(stream_id=stream_id, location=on_cells, &
      &          variable_name="cell_area", return_pointer=cell_data_broadcast)
    CALL read_2D(stream_id=stream_id, location=on_vertices, &
      &          variable_name="dual_area", return_pointer=vertex_data_broadcast)

    CALL closeFile(stream_id)

    CALL openinputfile(stream_id, dynamics_grid_filename(1), patch_2d, &
      &                read_netcdf_distribute_method)

    CALL read_2D(stream_id=stream_id, location=on_cells, &
      &          variable_name="cell_area", return_pointer=cell_data_distribute)
    CALL read_2D(stream_id=stream_id, location=on_vertices, &
      &          variable_name="dual_area", return_pointer=vertex_data_distribute)

    CALL closeFile(stream_id)

    IF ( MAXVAL(ABS(cell_data_distribute - cell_data_broadcast )) > 0.0_wp ) &
      CALL finish(method_name, "Cell data check failed")
    IF ( MAXVAL(ABS(vertex_data_distribute - vertex_data_broadcast )) > 0.0_wp ) &
      CALL finish(method_name, "Vertex data check failed")

    !---------------------------------------------------------------------

    RETURN 
    !---------------------------------------------------------------------

  END SUBROUTINE ocean_test_quads
  !-------------------------------------------------------------------------

END MODULE mo_ocean_testbed_quads

