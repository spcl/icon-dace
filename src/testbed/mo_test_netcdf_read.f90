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

! Tests the mo_read_interface methods

MODULE mo_test_netcdf_read

  USE mo_kind,                ONLY: wp
  USE mo_exception,           ONLY: message, finish
  USE mo_timer,               ONLY: ltimer, activate_sync_timers, timers_level
  USE mo_parallel_config,     ONLY: nproma
  USE mo_mpi,                 ONLY: my_process_is_mpi_workroot
  USE mo_communication,       ONLY: exchange_data

  USE mo_icon_testbed_config, ONLY: testfile_3D_time, testfile_2D_time

  USE mo_model_domain,        ONLY: t_patch, p_patch
  USE mo_atmo_model,          ONLY: construct_atmo_model, destruct_atmo_model
  USE mo_read_interface,      ONLY: openInputFile, closeFile, on_cells, &
    &                               t_stream_id, read_1D, read_2D_time, &
    &                               read_3D_time, read_netcdf_broadcast_method
  USE mo_netcdf_errhandler,   ONLY: nf
  USE mo_netcdf
!-------------------------------------------------------------------------
IMPLICIT NONE
PRIVATE

  PUBLIC :: test_netcdf_read, netcdf_open_output, netcdf_close_output
  PUBLIC :: netcdf_write_oncells_2D_time, netcdf_write_oncells_3D_time

  INTERFACE netcdf_write_oncells_2D_time
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_2D_time_filename
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_2D_time_fileid
  END INTERFACE netcdf_write_oncells_2D_time

  INTERFACE netcdf_write_oncells_3D_time
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_3D_time_filename
    MODULE PROCEDURE netcdf_write_REAL_ONCELLS_3D_time_fileid
  END INTERFACE netcdf_write_oncells_3D_time


CONTAINS
  
  !-------------------------------------------------------------------------
  !>
  !! Test reading amip aerosol
  !! tesetd with /pool/data/ICON/import/incoming/bc_aeropt_kinne/bc_aeropt_kinne_sw_b14_fin_2000.nc
  SUBROUTINE test_netcdf_read(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch), POINTER   :: patch
    REAL(wp), POINTER :: lnwl_array(:), levels_array(:), times_array(:)
    REAL(wp), POINTER :: aod(:,:,:,:), asy(:,:,:,:)   ! is (nproma, lnwl, blocks, time (months) )
    REAL(wp), POINTER :: z_aer_fine_mo(:,:,:,:)
    INTEGER :: levels, lnwl_size, file_id
    TYPE(t_stream_id) :: stream_id

    CHARACTER(*), PARAMETER :: method_name = "mo_test_netcdf_read:test_netcdf_read"

    !---------------------------------------------------------------------
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    patch => p_patch(1)
    levels  = patch%nlev
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! test sst
    CALL message(method_name,   testfile_3D_time(1))
    CALL openinputfile(file_id, filename = testfile_3D_time(1))

    !---------------------------------------------------------------------
    CALL read_1D(file_id=file_id, variable_name = "lnwl", &
      &          return_pointer=lnwl_array)
    CALL read_1D(file_id=file_id, variable_name = "lev", &
      &          return_pointer=levels_array)
    CALL read_1D(file_id=file_id, variable_name="time", &
      &          return_pointer=times_array)

    CALL closeFile(file_id)

    !---------------------------------------------------------------------
    CALL openinputfile(stream_id, testfile_3D_time(1), patch, &
      &                read_netcdf_broadcast_method)
    !-----------------------------------------------------
    ! example with non-allocated arrays
    ! mote that allocation time dim will be 1:12, as in the file
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name="aod", return_pointer=aod)

    !-----------------------------------------------------
    ! example with allocated arrays, time dim is 0:13
    lnwl_size = SIZE(lnwl_array, 1)
    ALLOCATE(asy(nproma, lnwl_size, patch%nblks_c, 0:13 ))
    ! read previous month, should be filled form the previous read though
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name="asy", fill_array=asy(:,:,:,0:0), &
      &               start_timestep=12, end_timestep=12, &
      &               levelsDimName = "lnwl") ! optional, just for checking
    ! read current year
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name="asy", fill_array=asy(:,:,:,1:12), &
      &               levelsDimName="lnwl") ! optional, just for checking
    ! read next month
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name="asy", fill_array=asy(:,:,:,13:13), &
      &               start_timestep=1, end_timestep=1)

    !-----------------------------------------------------
    ! check Z-axis vertical levels
    IF (SIZE(levels_array, 1) /= levels) THEN
      write(0,*) "patch levels=", levels, " SIZE(levels_array)=", SIZE(levels_array, 1)
      CALL finish(method_name, "SIZE(levels_array, 1) /= levels")
    ENDIF
    ALLOCATE(z_aer_fine_mo(nproma, levels, patch%nblks_c, 0:13 ))
    ! read previous month, should be filled form the previous read though
    ! read current year
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name="z_aer_fine_mo", &
      &               fill_array=z_aer_fine_mo(:,:,:,1:12), &
      &               levelsDimName = "lev") ! optional, just for checking

    !---------------------------------------------------------------------
    CALL closeFile(stream_id)
    DEALLOCATE(lnwl_array, levels_array, times_array, aod, asy, z_aer_fine_mo)

    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    CALL destruct_atmo_model()


  END SUBROUTINE test_netcdf_read
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  !!
  SUBROUTINE test_netcdf_read_demo_1(namelist_filename,shr_namelist_filename)

    CHARACTER(LEN=*), INTENT(in) :: namelist_filename
    CHARACTER(LEN=*), INTENT(in) :: shr_namelist_filename

    TYPE(t_patch), POINTER   :: patch
    REAL(wp), POINTER :: fill_3D_time_array(:,:,:,:)   ! is (nproma, vertical_levels, blocks, time (months))
    REAL(wp), POINTER :: fill_2D_time_array(:,:,:)  ! is (nproma, blocks, time )
    INTEGER :: return_status
    TYPE(t_stream_id) :: stream_id

    CHARACTER(*), PARAMETER :: method_name = "test_netcdf_read"

    !---------------------------------------------------------------------
    ltimer = .false.
    timers_level = 0
    activate_sync_timers = .false.
    CALL construct_atmo_model(namelist_filename,shr_namelist_filename)
    patch => p_patch(1)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! test sst
    CALL message(method_name,   testfile_2D_time(1))

    CALL openinputfile(stream_id, testfile_2D_time(1), patch, &
      &                read_netcdf_broadcast_method)

    !---------------------------------------------------------------------
    ! will read all timesteps in the file
    CALL read_2D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name=TRIM(testfile_2D_time(2)), &
      &               return_pointer=fill_2D_time_array)

    ! here the fill_2D_time_array is allocated and filled
    ! write the sst array to output ONLY for comparing
    return_status = netcdf_write_oncells_2d_time("testfile_2D_time_alltimes.nc", TRIM(testfile_2D_time(2)), &
      & fill_2D_time_array, patch)
    DEALLOCATE(fill_2D_time_array)
    !---------------------------------------------------------------------

    ! Example of reading selective timesteps with explicit shape
    ! say 13 months
    ALLOCATE(fill_2D_time_array(nproma, patch%nblks_c, 0:13 ))
    CALL read_2D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name=TRIM(testfile_2D_time(2)), &
      &               fill_array=fill_2D_time_array(:,:,0:), &
      &               start_timestep=12, end_timestep=12)
    CALL read_2D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name=TRIM(testfile_2D_time(2)), &
      &               fill_array=fill_2D_time_array(:,:,1:), &
      &               start_timestep=1,  end_timestep=12)
    CALL read_2D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name=TRIM(testfile_2D_time(2)), &
      &               fill_array=fill_2D_time_array(:,:,13:), &
      &               start_timestep=1, end_timestep=1)
    ! write(0,*) "shape of ", TRIM(testfile_2D_time(2)), ":", SHAPE(fill_2D_time_array)
    !---------------------------------------------------------------------
    ! write the sst array to output for comparing
    return_status = netcdf_write_oncells_2d_time("testfile_2D_14times.nc", TRIM(testfile_2D_time(2)), &
      & fill_2D_time_array, patch)
    CALL closeFile(stream_id)
    DEALLOCATE(fill_2D_time_array)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! test O3
    CALL message(method_name,   testfile_3D_time(1))
    CALL openinputfile(stream_id, testfile_3D_time(1), patch, &
      &                read_netcdf_broadcast_method)
    CALL read_3D_time(stream_id=stream_id, location=on_cells, &
      &               variable_name=TRIM(testfile_3D_time(2)), &
      &               return_pointer=fill_3D_time_array)
    !---------------------------------------------------------------------
    ! write the o3 array to output for comparing
    CALL message(method_name,   "write testfile_3D_time(1)")
    return_status = netcdf_write_oncells_3d_time(  &
      & filename=testfile_3D_time(1),              &
      & variable_name=TRIM(testfile_3D_time(2)),   &
      & write_array=fill_3D_time_array,            &
      & patch=patch)
    CALL closeFile(stream_id)
    DEALLOCATE(fill_3D_time_array)
    !---------------------------------------------------------------------
!    NULLIFY(fill_3D_time_array)
!    ALLOCATE(fill_3D_time_array(nproma, levels,  patch%nblks_c, 0:0 ))


    !---------------------------------------------------------------------
    ! Carry out the shared clean-up processes
    CALL destruct_atmo_model()
     

  END SUBROUTINE test_netcdf_read_demo_1

  !-------------------------------------------------------------------------

  INTEGER FUNCTION netcdf_open_output(filename)
    CHARACTER(LEN=*), INTENT(IN) :: filename

    INTEGER :: file_id, old_mode

    IF( my_process_is_mpi_workroot()  ) THEN
      CALL nf(nf90_create(TRIM(filename), &
        &                 IOR(nf90_clobber, nf90_64bit_offset), &
        &                 file_id), &
        &     TRIM(filename))
      CALL nf(nf90_set_fill(file_id, nf90_nofill, old_mode), TRIM(filename))
    ELSE
        file_id = -1 ! set it to an invalid value
    ENDIF

    netcdf_open_output = file_id

  END FUNCTION netcdf_open_output

  !-------------------------------------------------------------------------

  INTEGER FUNCTION netcdf_close_output(file_id)
    INTEGER, INTENT(IN) :: file_id

    netcdf_close_output = -1
    IF( my_process_is_mpi_workroot()  ) THEN
        netcdf_close_output = nf90_close(file_id)
    ENDIF

  END FUNCTION netcdf_close_output

  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_2D_time_filename(filename, variable_name, write_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: write_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'netcdf_write_REAL_ONCELLS_3D_time_filename'

    netcdf_write_REAL_ONCELLS_2D_time_filename = 0

    file_id = netcdf_open_output(filename)

    netcdf_write_REAL_ONCELLS_2D_time_filename = &
      & netcdf_write_REAL_ONCELLS_2D_time_fileid(file_id, variable_name, write_array, patch)

    return_status = netcdf_close_output(file_id)

  END FUNCTION netcdf_write_REAL_ONCELLS_2D_time_filename
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_3D_time_filename(filename, variable_name, write_array, patch)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp), POINTER            :: write_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch

    INTEGER :: file_id
    INTEGER :: return_status
    CHARACTER(LEN=*), PARAMETER :: method_name = 'netcdf_write_REAL_ONCELLS_3D_time_filename'

    netcdf_write_REAL_ONCELLS_3D_time_filename = 0

    file_id = netcdf_open_output(filename)

    netcdf_write_REAL_ONCELLS_3D_time_filename = &
      & netcdf_write_REAL_ONCELLS_3D_time_fileid(file_id, variable_name, write_array, patch)

    return_status = netcdf_close_output(file_id)

  END FUNCTION netcdf_write_REAL_ONCELLS_3D_time_filename
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  ! this function is meant only for checking the read methods
  ! Do not use for regular output
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_2D_time_fileid(file_id, variable_name, write_array, patch)
    INTEGER, INTENT(IN)  :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp),      POINTER       :: write_array(:,:,:)
    TYPE(t_patch), TARGET        :: patch

    REAL(wp), ALLOCATABLE        :: output_array(:,:)
    INTEGER :: total_number_of_cells, array_time_steps
    INTEGER :: dim_number_of_cells, dim_array_time_steps
    INTEGER :: output_shape(2), dim_write_shape(2), diff_shape(2)
    INTEGER :: variable_id
    INTEGER :: start_timestep, end_timestep, timestep

!    INTEGER                      :: i,j,t
    CHARACTER(LEN=*), PARAMETER  :: method_name = 'netcdf_write_REAL_ONCELLS_3D_time_fileid'

    netcdf_write_REAL_ONCELLS_2D_time_fileid = 0

    start_timestep = LBOUND(write_array, 3)
    end_timestep   = UBOUND(write_array, 3)

    ALLOCATE(output_array(MERGE(patch%n_patch_cells_g, 0, &
      &                         my_process_is_mpi_workroot()), &
      &                   start_timestep:end_timestep))
    DO timestep = start_timestep, end_timestep
      CALL exchange_data(in_array=write_array(:,:,timestep), &
        &                out_array=output_array(:,timestep), &
        &                gather_pattern=patch%comm_pat_gather_c)
    END DO

    !----------------------------------------------------------------------
    ! Write only from mpi_workroot
    IF( my_process_is_mpi_workroot()  ) THEN

      ! Write Dimensions
      total_number_of_cells = patch%n_patch_cells_g
      array_time_steps      = SIZE(write_array, 3)
      output_shape          = (/ total_number_of_cells, array_time_steps /)
      diff_shape            = (SHAPE(output_array) - output_shape)

      IF ( MAXVAL(ABS( diff_shape)) /= 0 ) THEN
        WRITE(0,*) " gather array shape:", SHAPE(output_array)
        WRITE(0,*) " computed array shape:", output_shape
        CALL finish(method_name, "gather array shape is nor correct")
      ENDIF

      CALL nf(nf90_def_dim(file_id, 'ncells',  total_number_of_cells,  dim_number_of_cells),  variable_name)
      CALL nf(nf90_def_dim(file_id, 'time',    array_time_steps,       dim_array_time_steps), variable_name)
      dim_write_shape = (/ dim_number_of_cells, dim_array_time_steps /)

      ! define variable
      ! WRITE(0,*) " define variable:", variable_name
!      CALL nf(nf90_def_var(file_id, variable_name, nf90_double, dim_write_shape, &
!        & variable_id), variable_name)
      CALL nf(nf90_def_var(file_id, variable_name, nf90_float, dim_write_shape, &
        & variable_id), variable_name)

      CALL nf(nf90_enddef(file_id), variable_name)

      ! WRITE(0,*) " write array..."
      CALL nf(nf90_put_var(file_id, variable_id, output_array(:,:)), variable_name)
      ! CALL nf(nf90_put_var(file_id, variable_id, REAL(output_array(:,:,:))), variable_name)

    ENDIF

    DEALLOCATE(output_array)

  END FUNCTION netcdf_write_REAL_ONCELLS_2D_time_fileid
  !-------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  !>
  ! this function is meant only for checking the read methods
  ! Do not use for regular output
  INTEGER FUNCTION netcdf_write_REAL_ONCELLS_3D_time_fileid(file_id, variable_name, write_array, patch)
    INTEGER, INTENT(IN)  :: file_id
    CHARACTER(LEN=*), INTENT(IN) :: variable_name
    REAL(wp),      POINTER       :: write_array(:,:,:,:)
    TYPE(t_patch), TARGET        :: patch

    REAL(wp), ALLOCATABLE        :: output_array(:,:,:)
    INTEGER :: total_number_of_cells, array_vertical_levels, array_time_steps
    INTEGER :: dim_number_of_cells, dim_vertical_levels, dim_array_time_steps
    INTEGER :: output_shape(3), dim_write_shape(3), diff_shape(3)
    INTEGER :: variable_id
    INTEGER :: start_timestep, end_timestep, timestep, nlev

!    INTEGER                      :: i,j,t
    CHARACTER(LEN=*), PARAMETER  :: method_name = 'netcdf_write_REAL_ONCELLS_3D_time_fileid'

    netcdf_write_REAL_ONCELLS_3D_time_fileid = 0

    nlev           = SIZE(write_array, 2)
    start_timestep = LBOUND(write_array, 4)
    end_timestep   = UBOUND(write_array, 4)

    ALLOCATE(output_array(MERGE(patch%n_patch_cells_g, 0, &
      &                         my_process_is_mpi_workroot()), &
      &                   nlev, start_timestep:end_timestep))
    DO timestep = start_timestep, end_timestep
      CALL exchange_data(in_array=write_array(:,:,:, timestep), &
        &                out_array=output_array(:, :, timestep), &
        &                gather_pattern=patch%comm_pat_gather_c)
    END DO

    !----------------------------------------------------------------------
    ! Write only from mpi_workroot
    IF( my_process_is_mpi_workroot()  ) THEN

      ! Write Dimensions
      total_number_of_cells = patch%n_patch_cells_g
      array_vertical_levels = SIZE(write_array, 2)
      array_time_steps      = SIZE(write_array, 4)
      output_shape          = (/ total_number_of_cells, array_vertical_levels, array_time_steps /)
      diff_shape            = (SHAPE(output_array) - output_shape)

      IF ( MAXVAL(ABS( diff_shape)) /= 0 ) THEN
        WRITE(0,*) " gather array shape:", SHAPE(output_array)
        WRITE(0,*) " computed array shape:", output_shape
        CALL finish(method_name, "gather array shape is nor correct")
      ENDIF

      CALL nf(nf90_def_dim(file_id, 'ncells',  total_number_of_cells,  dim_number_of_cells),  variable_name)
      CALL nf(nf90_def_dim(file_id, 'plev',    array_vertical_levels,  dim_vertical_levels),  variable_name)
      CALL nf(nf90_def_dim(file_id, 'time',    array_time_steps,       dim_array_time_steps), variable_name)
      dim_write_shape = (/ dim_number_of_cells, dim_vertical_levels,  dim_array_time_steps /)

      ! define variable
      ! WRITE(0,*) " define variable:", variable_name
!      CALL nf(nf90_def_var(file_id, variable_name, nf90_double, dim_write_shape, &
!        & variable_id), variable_name)
      CALL nf(nf90_def_var(file_id, variable_name, nf90_float, dim_write_shape, &
        & variable_id), variable_name)

      CALL nf(nf90_enddef(file_id), variable_name)

!      DO t=1, array_time_steps
!        DO i=1, total_number_of_cells
!          DO j=1, array_vertical_levels
!            write(0,*) t,i,j, ":", output_array(i,j,t)
!          ENDDO
!        ENDDO
!      ENDDO
      ! write array
      ! WRITE(0,*) " write array..."
      CALL nf(nf90_put_var(file_id, variable_id, output_array(:,:,:)), variable_name)
      ! CALL nf(nf90_put_var(file_id, variable_id, REAL(output_array(:,:,:))), variable_name)

    ENDIF

    DEALLOCATE(output_array)

  END FUNCTION netcdf_write_REAL_ONCELLS_3D_time_fileid
  !-------------------------------------------------------------------------


END MODULE mo_test_netcdf_read

