! Contains namelists for parallel run control.
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

MODULE mo_parallel_nml

  USE mo_io_units,            ONLY: nnml, nnml_output
  USE mo_exception,           ONLY: finish, message, message_text
  USE mo_namelist,            ONLY: position_nml, POSITIONED, open_nml, close_nml
  USE mo_master_control,      ONLY: use_restart_namelists
  USE mo_mpi,                 ONLY: my_process_is_stdio
  USE mo_restart_nml_and_att, ONLY: open_tmpfile, store_and_close_namelist, open_and_restore_namelist, close_tmpfile
  USE mo_io_units,            ONLY: filename_max
  USE mo_impl_constants,      ONLY: max_dom
  USE mo_nml_annotate,        ONLY: temp_defaults, temp_settings

  USE mo_parallel_config,     ONLY: &
    & config_n_ghost_rows        => n_ghost_rows,        &
    & config_division_method     => division_method,     &
    & config_division_file_name  => division_file_name,  &
    & config_write_div_to_file   => write_div_to_file,   &
    & config_use_div_from_file   => use_div_from_file,   &
    & config_ldiv_phys_dom       => ldiv_phys_dom,       &
    & config_l_log_checks        => l_log_checks,        &
    & config_l_fast_sum          => l_fast_sum,          &
    & config_p_test_run          => p_test_run,          &
    & config_num_test_pe         => num_test_pe,         &
    & config_l_test_openmp       => l_test_openmp,       &
    & config_num_restart_procs   => num_restart_procs,   &
    & config_num_io_procs        => num_io_procs,        &
    & config_num_io_procs_radar        => num_io_procs_radar,        &
    & config_pio_type            => pio_type,            &
    & config_num_prefetch_proc   => num_prefetch_proc,   &
    & config_proc0_shift         => proc0_shift,         &
    & config_iorder_sendrecv     => iorder_sendrecv,     &
    & set_nproma, &
    & set_nproma_nblocks, &
    & set_nproma_nblocks_sub, &
    & config_nblocks_c           => nblocks_c,             &
    & config_nblocks_e           => nblocks_e,             &
    & config_use_nblocks_c       => ignore_nproma_use_nblocks_c,     &
    & config_use_nblocks_e       => ignore_nproma_use_nblocks_e,     &
    & config_nproma_sub          => nproma_sub,            &
    & config_nblocks_sub         => nblocks_sub,           &
    & config_use_nblocks_sub     => ignore_nproma_sub_use_nblocks_sub,     &
    & config_use_icon_comm       => use_icon_comm,        &
    & config_icon_comm_debug     => icon_comm_debug,        &
    & div_geometric, check_parallel_configuration,          &
    & config_max_sr_buffer_size => max_send_recv_buffer_size, &
    & config_use_omp_input          => use_omp_input,         &
    & config_use_dycore_barrier => use_dycore_barrier,        &
    & config_itype_exch_barrier => itype_exch_barrier,        &
    & config_use_dp_mpi2io      => use_dp_mpi2io,             &
    & config_icon_comm_method   => icon_comm_method,          &
    & config_max_no_of_comm_var => max_no_of_comm_variables,  &
    & config_max_no_of_comm_proc => max_no_of_comm_processes, &
    & config_max_no_of_comm_patt => max_no_of_comm_patterns,  &
    & config_sync_barrier_mode   => sync_barrier_mode,        &
    & config_max_mpi_message_size => max_mpi_message_size,    &
    & config_use_physics_barrier  => use_physics_barrier,     &
    & config_restart_chunk_size => restart_chunk_size,        &
    & config_restart_load_scale_max => restart_load_scale_max, &
    & config_io_proc_chunk_size => io_proc_chunk_size,        &
    & config_num_dist_array_replicas => num_dist_array_replicas, &
    & config_io_process_stride => io_process_stride,          &
    & config_io_process_rotate => io_process_rotate,          &
    & comm_pattern_type_orig,                                 &
    & config_default_comm_pattern_type => default_comm_pattern_type

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: read_parallel_namelist


  CONTAINS

  !-------------------------------------------------------------------------
  !
  SUBROUTINE read_parallel_namelist( filename )

    ! ------------------------------------------------------------------------
    ! Number of rows of ghost cells
    INTEGER :: n_ghost_rows

    INTEGER :: division_method(0:max_dom)
!                      div_geometric = 1  ! Geometric subdivision
!                      ext_div_from_file = 201 ! Read from file

    CHARACTER(LEN=filename_max) :: division_file_name(0:max_dom) ! if ext_div_from_file

    LOGICAL :: write_div_to_file
    LOGICAL :: use_div_from_file

    ! Flag if (in case of merged domains) physical domains shall be considered for
    ! computing the domain decomposition
    LOGICAL :: ldiv_phys_dom

    ! Flag if checks in a verification run should be logged
    LOGICAL :: l_log_checks

    ! Flag if fast but nonreproducible sum should be used
    LOGICAL :: l_fast_sum

    ! Please note for the following variables: The default settings are for NO_MPI runs!

    ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
    ! model whereas the other PEs do a real parallelized run
    LOGICAL :: p_test_run

    ! use more than 1 PE for verification if p_test_run and num_test_pe is set
    ! to a value > 1
    INTEGER :: num_test_pe
    LOGICAL :: use_dycore_barrier ! put an mpi barrier before the dycore to synchronize MPI tasks
    LOGICAL :: use_physics_barrier
    INTEGER :: itype_exch_barrier ! 1: put an mpi barrier at the beginning of exchange calls to synchronize MPI tasks
                                  ! 2: put an mpi barrier after MPI_WAIT to synchronize MPI tasks
                                  ! 3: 1+2

    ! if l_test_openmp is set together with p_test_run, then the verification PE uses
    ! only 1 thread. This allows for verifying the OpenMP implementation
    LOGICAL :: l_test_openmp

    LOGICAL :: use_icon_comm
    LOGICAL :: icon_comm_debug
    INTEGER :: max_send_recv_buffer_size
    INTEGER :: max_mpi_message_size
    INTEGER :: icon_comm_method
    INTEGER :: max_no_of_comm_variables
    INTEGER :: max_no_of_comm_processes
    INTEGER :: max_no_of_comm_patterns
    INTEGER :: sync_barrier_mode

    ! Type of parallel I/O
    INTEGER :: pio_type
    INTEGER :: num_io_procs
    INTEGER :: num_io_procs_radar

    ! Number of restart PEs (0 means, worker 0 writes restart (to be backward compatible)
    INTEGER :: num_restart_procs

    ! Number of PEs used for async prefetching of input (0 means, the worker PE0 prefetches lateral boundary input)
    INTEGER :: num_prefetch_proc

    ! Shift of processor 0 in domain decomposition, e.g. to use proc 0 for input only
    INTEGER :: proc0_shift

    ! Use OpenMP-parallelized input for atmospheric input data (in initicon), 
    ! i.e. overlapping of reading data, communicating data and computing statistics
    LOGICAL :: use_omp_input

    ! Order of send/receive sequence in exchange routines
    ! 1 = irecv, send
    ! 2 = isend, recv
    ! 3 = irecv, isend
    INTEGER :: iorder_sendrecv

    INTEGER :: nproma    ! inner loop length/vector length
    INTEGER :: nblocks_c   ! inner loop number of blocks used for cells
    INTEGER :: nblocks_e   ! inner loop number of blocks used for edges

    INTEGER :: nproma_sub    ! secondary nproma, by default the same as the basic one
    INTEGER :: nblocks_sub   ! number of blocks to subdivide the basic nproma into

    LOGICAL :: use_dp_mpi2io

    ! The (asynchronous) restart is capable of writing and communicating
    ! more than one 2D slice at once
    INTEGER :: restart_chunk_size
    ! The multifile checkpointing framework is capable of reading and distributing
    ! full 3d arrays (if there are less than restart_load_scale_max work PE per file.
    INTEGER :: restart_load_scale_max

    ! The (asynchronous) name list output is capable of writing and communicating
    ! more than one 2D slice at once
    INTEGER :: io_proc_chunk_size

    INTEGER :: io_process_stride, io_process_rotate

    ! number of replications being stored in the distributed arrays of the
    ! t_patch_pre
    INTEGER :: num_dist_array_replicas

    ! default implementation of mo_communication to be used
    ! 1 = comm_pattern_type_orig
    ! 2 = comm_pattern_type_yaxt
    INTEGER :: default_comm_pattern_type

    NAMELIST /parallel_nml/ n_ghost_rows,  division_method, ldiv_phys_dom, &
      & l_log_checks,      l_fast_sum,          &
      & p_test_run, num_test_pe, l_test_openmp,       &
      & num_restart_procs, proc0_shift,         &
      & num_io_procs,      pio_type,            &
      & num_io_procs_radar,                     &
      & iorder_sendrecv,                        &
      & nproma, nblocks_c, nblocks_e,           &
      & nproma_sub, nblocks_sub,                &
      & use_icon_comm, &
      & icon_comm_debug, max_send_recv_buffer_size, &
      & division_file_name, use_dycore_barrier, &
      & write_div_to_file, use_div_from_file, &
      & use_dp_mpi2io, itype_exch_barrier,                &
      & icon_comm_method, max_no_of_comm_variables,       &
      & max_no_of_comm_processes, max_no_of_comm_patterns, &
      & sync_barrier_mode, max_mpi_message_size, use_physics_barrier, &
      & restart_chunk_size, io_proc_chunk_size, num_prefetch_proc, &
      & num_dist_array_replicas, io_process_stride, io_process_rotate, &
      & default_comm_pattern_type, use_omp_input, restart_load_scale_max

    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER :: istat
    INTEGER :: funit
    INTEGER :: iunit
    CHARACTER(len=*), PARAMETER ::   &
            &  routine = 'mo_parallel_nml:read_parallel_nml'

    !--------------------------------------------
    ! set default values
    !--------------------------------------------
    ! Number of rows of ghost cells
    n_ghost_rows = 1
    division_method(:) = div_geometric
    division_file_name(:) = ""
    write_div_to_file = .FALSE.
    use_div_from_file = .FALSE.

    ! Flag if (in case of merged domains) physical domains shall be considered for
    ! computing the domain decomposition
    ldiv_phys_dom = .TRUE.

    ! Flag if checks in a verification run should be logged
    l_log_checks = .FALSE.

    ! Flag if fast but nonreproducible sum should be used
    l_fast_sum = .FALSE.

    ! Please note for the following variables: The default settings are for NO_MPI runs!
    ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
    ! model whereas the other PEs do a real parallelized run
    p_test_run = .FALSE.
    num_test_pe = -1

    ! The barriers should be used for dedicated tests only, not for production runs
    use_dycore_barrier = config_use_dycore_barrier
    use_physics_barrier= config_use_physics_barrier
    itype_exch_barrier = 0

    ! if l_test_openmp is set together with p_test_run, then the verification PE uses
    ! only 1 thread. This allows for verifying the OpenMP implementation
    l_test_openmp = .FALSE.

    use_icon_comm = .FALSE.
    icon_comm_debug = .FALSE.
    max_send_recv_buffer_size = config_max_sr_buffer_size
    max_mpi_message_size      = config_max_mpi_message_size
    icon_comm_method          = config_icon_comm_method
    max_no_of_comm_variables  = config_max_no_of_comm_var
    max_no_of_comm_processes  = config_max_no_of_comm_proc
    max_no_of_comm_patterns   = config_max_no_of_comm_patt
    sync_barrier_mode = 0

    ! Type of parallel I/O
    pio_type = 1
    num_io_procs = 0
    num_io_procs_radar = 0

    ! Number of restart output PEs; if 0, worker 0 does the work
    num_restart_procs = 0

    ! The number of PEs used for async prefetching of input (0 means, the worker PE0 prefetches input)
    num_prefetch_proc = 1

    ! Shift of processor 0 in domain decomposition, set to 1 in order to use proc 0 for input only
    proc0_shift = 0

    ! Use OpenMP-parallelized input for atmospheric input data (in initicon), 
    ! i.e. overlapping of reading data, communicating data and computing statistics
    use_omp_input = .FALSE.

    ! Order of send/receive sequence in exchange routines
    ! 1 = irecv, send
    ! 2 = isend, recv
    ! 3 = irecv, isend
    iorder_sendrecv = 1

    ! inner loop length/vector length
    nproma = 0
    nblocks_c = 0
    nblocks_e = 0

    ! secondary inner loop length/vector length
    nproma_sub = 0
    nblocks_sub = 0

    ! MPI gather to output processes in DOUBLE PRECISION
    use_dp_mpi2io = .FALSE.

    restart_chunk_size = -1
    restart_load_scale_max = 1

    io_proc_chunk_size = -1

    num_dist_array_replicas = 1

    io_process_stride = -1
    io_process_rotate = 0
    default_comm_pattern_type = comm_pattern_type_orig

    !----------------------------------------------------------------
    ! If this is a resumed integration, overwrite the defaults above
    ! by values in the previous integration.
    !----------------------------------------------------------------
    IF (use_restart_namelists()) THEN
      funit = open_and_restore_namelist('parallel_nml')
      READ(funit,NML=parallel_nml)
      CALL close_tmpfile(funit)
    END IF

    !--------------------------------------------------------------------
    ! Read user's (new) specifications (Done so far by all MPI processes)
    !--------------------------------------------------------------------
    CALL open_nml(TRIM(filename))
    CALL position_nml ('parallel_nml', STATUS=istat)
    IF (my_process_is_stdio()) THEN
      iunit = temp_defaults()
      WRITE(iunit, parallel_nml)     ! write defaults to temporary text file
    END IF
    SELECT CASE (istat)
    CASE (POSITIONED)
      READ (nnml, parallel_nml)                                         ! overwrite default settings
      IF (my_process_is_stdio()) THEN
        iunit = temp_settings()
        WRITE(iunit, parallel_nml)     ! write settings to temporary text file
      END IF
    END SELECT
    CALL close_nml


    !-----------------------------------------------------
    ! Store the namelist for restart
    !-----------------------------------------------------
    IF(my_process_is_stdio())  THEN
      funit = open_tmpfile()
      WRITE(funit,NML=parallel_nml)
      CALL store_and_close_namelist(funit, 'parallel_nml')
    ENDIF
    !-----------------------------------------------------
    ! Write the final namelist to an ASCII file
    IF (my_process_is_stdio()) WRITE(nnml_output,nml=parallel_nml)

    !-----------------------------------------------------
    ! fill_parallel_nml_configure
    config_n_ghost_rows        = n_ghost_rows
    config_division_method(:)  = division_method(:)
    config_division_file_name(:) = division_file_name(:)
    config_write_div_to_file   = write_div_to_file
    config_use_div_from_file   = use_div_from_file
    config_ldiv_phys_dom       = ldiv_phys_dom
    config_l_log_checks        = l_log_checks
    config_l_fast_sum          = l_fast_sum
    config_p_test_run          = p_test_run
    config_num_test_pe         = num_test_pe
    config_l_test_openmp       = l_test_openmp
    config_num_restart_procs   = num_restart_procs
    config_num_io_procs        = num_io_procs
    config_num_io_procs_radar  = num_io_procs_radar
    config_pio_type            = pio_type
    config_num_prefetch_proc   = num_prefetch_proc
    config_proc0_shift         = proc0_shift
    config_use_omp_input       = use_omp_input
    config_iorder_sendrecv     = iorder_sendrecv
    CALL set_nproma_nblocks(nproma, nblocks_c, nblocks_e)
    CALL set_nproma_nblocks_sub(nproma_sub, nblocks_sub)

    config_use_icon_comm       = use_icon_comm
    config_icon_comm_debug     = icon_comm_debug
    config_icon_comm_method    = icon_comm_method
    config_max_no_of_comm_var  = max_no_of_comm_variables
    config_max_no_of_comm_proc = max_no_of_comm_processes
    config_max_no_of_comm_patt = max_no_of_comm_patterns
    config_sync_barrier_mode   = sync_barrier_mode
    config_max_sr_buffer_size   = max_send_recv_buffer_size
    config_max_mpi_message_size = max_mpi_message_size
    config_use_dycore_barrier   = use_dycore_barrier
    config_use_physics_barrier  = use_physics_barrier
    config_itype_exch_barrier   = itype_exch_barrier
    config_use_dp_mpi2io        = use_dp_mpi2io
    config_restart_chunk_size   = restart_chunk_size
    config_io_proc_chunk_size   = io_proc_chunk_size
    config_restart_load_scale_max = restart_load_scale_max
    config_num_dist_array_replicas = num_dist_array_replicas
    config_io_process_stride    = io_process_stride
    config_io_process_rotate    = io_process_rotate
    config_default_comm_pattern_type = default_comm_pattern_type
    !-----------------------------------------------------
    CALL check_parallel_configuration()

  END SUBROUTINE read_parallel_namelist
  !-------------------------------------------------------------------------

END MODULE mo_parallel_nml      
