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

MODULE mo_parallel_config

  USE mo_exception,          ONLY: message, finish, warning
  USE mo_io_units,           ONLY: filename_max
  USE mo_impl_constants,     ONLY: max_dom, max_num_io_procs, pio_type_async
  USE mo_util_string,        ONLY: int2string
  USE mo_mpi,                ONLY: p_bcast, process_mpi_all_workroot_id, &
  &                                process_mpi_all_comm, my_process_is_io

  IMPLICIT NONE

  PRIVATE
  ! Exported variables:
  PUBLIC :: nproma, nblocks_c, ignore_nproma_use_nblocks_c
  PUBLIC :: nproma_sub, nblocks_sub, ignore_nproma_sub_use_nblocks_sub
!
  PUBLIC :: n_ghost_rows,                                     &
       &  div_geometric, division_method, division_file_name,       &
       &  l_log_checks, l_fast_sum,   &
       &  ldiv_phys_dom, p_test_run, num_test_pe, l_test_openmp,    &
       &  pio_type, iorder_sendrecv, num_io_procs,                  &
       &  num_restart_procs, num_prefetch_proc, num_io_procs_radar, &
       &  use_icon_comm, icon_comm_debug, max_send_recv_buffer_size,&
       &  use_dycore_barrier, itype_exch_barrier, use_dp_mpi2io,    &
       &  icon_comm_method, icon_comm_openmp, max_no_of_comm_variables, &
       &  max_no_of_comm_processes, max_no_of_comm_patterns,        &
       &  sync_barrier_mode, max_mpi_message_size, use_physics_barrier, &
       &  restart_chunk_size, restart_load_scale_max, ext_div_from_file, &
       &  write_div_to_file, use_div_from_file, io_proc_chunk_size, &
       &  num_dist_array_replicas, comm_pattern_type_orig,          &
       &  comm_pattern_type_yaxt, default_comm_pattern_type,        &
       &  io_process_stride, io_process_rotate, proc0_shift,        &
       &  use_omp_input

  PUBLIC :: set_nproma, get_nproma, cpu_min_nproma, proc0_offloading, &
       &    check_parallel_configuration, use_async_restart_output, blk_no, idx_no, idx_1d,    &
       &    update_nproma_for_io_procs

  ! computing setup
  ! ---------------
  INTEGER  :: nproma = 1              ! inner loop length/vector length
  !$ACC DECLARE COPYIN(nproma)
  INTEGER  :: nblocks_c = 0
  LOGICAL  :: ignore_nproma_use_nblocks_c = .FALSE.

  ! Secondary nproma for, e.g., radiation chunking
  INTEGER  :: nproma_sub = -1
  !$ACC DECLARE COPYIN(nproma_sub)
  INTEGER  :: nblocks_sub = -1
  LOGICAL  :: ignore_nproma_sub_use_nblocks_sub = .TRUE.

  ! Number of rows of ghost cells
  INTEGER :: n_ghost_rows = 1

  ! Division method for area subdivision
  INTEGER, PARAMETER :: div_geometric = 1  ! Geometric subdivision
  INTEGER, PARAMETER :: ext_div_from_file = 201 ! Read from file

  INTEGER :: division_method(0:max_dom) = 1
  CHARACTER(LEN=filename_max) :: division_file_name(0:max_dom)! if ext_div_from_file
  LOGICAL :: use_div_from_file = .FALSE. ! check for domain decomposition from file
                                         ! if file is not available use division_method
                                         ! to generate decomposition online
  LOGICAL :: write_div_to_file = .FALSE. ! write result of domain decomposition to file

  ! Flag if (in case of merged domains) physical domains shall be considered for
  ! computing the domain decomposition
  LOGICAL :: ldiv_phys_dom = .FALSE.

  ! Flag if checks in a verification run should be logged
  LOGICAL :: l_log_checks = .false.

  ! Flag if fast but nonreproducible sum should be used
  LOGICAL :: l_fast_sum = .false.

  ! Please note for the following variables: The default settings are for NO_MPI runs!

  ! p_test_run indicates a verification run, i.e. a run where 1 PE runs the complete
  ! model whereas the other PEs do a real parallelized run
  LOGICAL :: p_test_run = .false.

  ! use more than 1 PE for verification if p_test_run and num_test_pe is set
  ! to a value > 1
  INTEGER :: num_test_pe

  LOGICAL :: use_dycore_barrier = .false. ! acivate an mpi barrier before the dycore
                                          ! to synchronize MPI tasks
  LOGICAL :: use_physics_barrier = .false. ! activate mpi barrier after the physics

  INTEGER :: itype_exch_barrier = 0  ! 1: put an mpi barrier at the beginning of exchange calls to synchronize MPI tasks
                                     ! 2: put an mpi barrier after MPI_WAIT to synchronize MPI tasks
                                     ! 3: 1+2

  ! if l_test_openmp is set together with p_test_run, then the verification PE uses
  ! only 1 thread. This allows for verifying the OpenMP implementation
  LOGICAL :: l_test_openmp = .false.

  LOGICAL :: use_icon_comm = .false.
  LOGICAL :: icon_comm_debug= .false.
  INTEGER :: max_send_recv_buffer_size = 262144  ! size in doubles (x8)
  INTEGER :: max_mpi_message_size      = 65536   ! size in doubles (x8)
  INTEGER :: max_no_of_comm_variables  = 64
  INTEGER :: max_no_of_comm_processes  = 64
  INTEGER :: max_no_of_comm_patterns   = 32
  INTEGER :: icon_comm_method = 1
  INTEGER :: sync_barrier_mode = 0

  LOGICAL :: icon_comm_openmp = .false.

  ! Flag whether async restart output is used, it is set in the main program:
  LOGICAL :: use_async_restart_output = .FALSE.

  ! Type of parallel I/O
  INTEGER :: pio_type = pio_type_async

  INTEGER :: num_io_procs = 0

  ! Output procs for radar forward operator
  INTEGER :: num_io_procs_radar = 0

  ! The number of PEs used for writing restart files (0 means, the worker PE0 writes)
  INTEGER :: num_restart_procs = 0

  ! The number of PEs used for async prefetching of input (0 means, the PE0 prefetches input)
  INTEGER :: num_prefetch_proc = 0

  ! Shift of processor 0 in domain decomposition, e.g. to use proc 0 for input only
  INTEGER :: proc0_shift = 0

  ! Derived variable to indicate hybrid mode with proc 0 doing stdio only
  LOGICAL :: proc0_offloading

  ! Use OpenMP-parallelized input for atmospheric input data (in initicon), 
  ! i.e. overlapping of reading data, communicating data and computing statistics
  LOGICAL :: use_omp_input = .FALSE.

  ! Order of send/receive sequence in exchange routines
  ! 1 = irecv, send
  ! 2 = isend, recv
  ! 3 = irecv, isend
  INTEGER :: iorder_sendrecv = 1

  ! Flag. Enable this flag if output fields shall be gathered in
  ! DOUBLE PRECISION. The resulting files are identical to the "normal"
  ! operation with simple REAL, since NetCDFs are written in single precision
  ! anyway by default and GRIB2 has 16 or 24 bit per data word.
  !
  LOGICAL :: use_dp_mpi2io

  ! The asynchronous and multifile checkpointing frameworks are capable of writing
  ! and communicating more than one 2D slice at once.
  INTEGER :: restart_chunk_size

  ! The multifile checkpointing framework is capable of reading and distributing
  ! full 3d arrays (if there are less than restart_load_scale_max work PE per file.
  INTEGER :: restart_load_scale_max = 1

  ! The (asynchronous) name list output is capable of writing and communicating
  ! more than one 2D slice at once
  INTEGER :: io_proc_chunk_size

  ! number of replications being stored in the distributed arrays of the
  ! t_patch_pre
  INTEGER :: num_dist_array_replicas

  ! use every nth process to do distributed netcdf reads
  INTEGER :: io_process_stride

  ! shift ranks doing I/O by this number
  INTEGER :: io_process_rotate

  ! switch between different implementations of mo_communication
  INTEGER, PARAMETER :: comm_pattern_type_orig = 1
  INTEGER, PARAMETER :: comm_pattern_type_yaxt = 2
  INTEGER :: default_comm_pattern_type

CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE check_parallel_configuration()

!   !local variables
!     INTEGER :: i_status, my_color, peer_comm, p_error
    CHARACTER(*), PARAMETER :: method_name = "check_parallel_configuration"

    !------------------------------------------------------------
    !  check the consistency of the parameters
    !------------------------------------------------------------
    IF (ignore_nproma_use_nblocks_c) THEN
      ! Sanity Check. (nblocks_c<=0) should never be the case here.
      IF (nblocks_c<=0) CALL finish(TRIM(method_name),'"nblocks_c" must be positive')
    ELSE
      IF (nproma<=0)  CALL finish(TRIM(method_name),'"nproma" must be positive')
    ENDIF
#if !defined (__SX__) && !defined (__NEC_VH__) && !defined(_OPENACC)
    ! migration helper: catch nproma's that were obviously intended
    !                   for a vector or GPU machine.
    IF (nproma>256) CALL warning(TRIM(method_name),'The value of "nproma" seems to be set for a vector machine!')
#endif

    IF (proc0_shift < 0) CALL finish(TRIM(method_name),'"proc0_shift" currently must be 0 (disable) or a positive number')
    proc0_offloading = proc0_shift > 0

    icon_comm_openmp = .false.
! check l_test_openmp
#ifndef _OPENMP
    IF (l_test_openmp) THEN
      CALL message(method_name, &
         & 'l_test_openmp has no effect if the model is compiled without OpenMP support')
      CALL message(method_name, &
         & '--> l_test_openmp set to .FALSE.')
      l_test_openmp = .FALSE.
    END IF
#else
    IF (icon_comm_method > 100) &
      & icon_comm_openmp = .true.
#endif

    ! check p_test_run, num_io_procs, num_restart_procs and num_prefetch_proc
#ifdef NOMPI
    ! Unconditionally set p_test_run to .FALSE. and num_io_procs to 0, num_restart_procs to 0
    ! and num_prefetch_proc to 0
    ! all other variables are already set correctly
    IF (p_test_run) THEN
      CALL message(method_name, &
       & 'p_test_run has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> p_test_run set to .FALSE.')
      p_test_run = .FALSE.
    END IF
    IF (num_io_procs /= 0) THEN
      CALL message(method_name, &
       & 'num_io_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> num_io_procs set to 0')
      num_io_procs = 0
    END IF
    IF (num_restart_procs /= 0) THEN
      CALL message(method_name, &
       & 'num_restart_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> num_restart_procs set to 0')
      num_restart_procs = 0
    END IF
    IF (num_prefetch_proc /= 0) THEN
      CALL message(method_name, &
       & 'num_prefetch_proc has no effect if the model is compiled with the NOMPI compiler directive')
      CALL message(method_name, &
       & '--> num_prefetch_proc set to 0')
      num_prefetch_proc = 0
    END IF

#else

    ! check n_ghost_rows
    IF (n_ghost_rows<1) THEN
      CALL finish(method_name, &
          & 'n_ghost_rows<1 in parallel_nml namelist is not allowed')
    END IF

    ! for safety only
    IF (num_io_procs < 0)      num_io_procs = 0
    IF (num_restart_procs < 0) num_restart_procs = 0
    IF (num_io_procs > MAX_NUM_IO_PROCS) THEN
      CALL finish(method_name, "Namelist parameter num_io_procs chosen too large ( > "//TRIM(int2string(MAX_NUM_IO_PROCS))//")!")
    END IF
    IF(num_prefetch_proc < 0) num_prefetch_proc = 0
    IF(num_prefetch_proc > 1) &
         CALL finish(TRIM(method_name),'The no of prefetch processor can be zero or one, but should not be set more than one!')

#endif

  END SUBROUTINE check_parallel_configuration
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE set_nproma(new_nproma)
    INTEGER, INTENT(IN) :: new_nproma

    nproma = new_nproma
    !$ACC UPDATE DEVICE(nproma) ASYNC(1) IF_PRESENT

  END SUBROUTINE set_nproma
  !-------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !>
  SUBROUTINE update_nproma_for_io_procs(nproma_max)
    INTEGER, INTENT(IN) :: nproma_max

    !adapt nproma for I/O procs
    CALL p_bcast(nproma_max, process_mpi_all_workroot_id, process_mpi_all_comm)
    IF (my_process_is_io()) THEN
      CALL set_nproma(nproma_max)
    ENDIF
  
  END SUBROUTINE update_nproma_for_io_procs
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  INTEGER FUNCTION get_nproma()

    get_nproma = nproma

  END FUNCTION get_nproma
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  FUNCTION cpu_min_nproma(nproma, min_nproma) RESULT(new_nproma)
    INTEGER, INTENT(IN) :: nproma, min_nproma
    INTEGER             :: new_nproma

#ifdef _OPENACC
    new_nproma = nproma
#else
    new_nproma = MIN(nproma, min_nproma)
#endif

  END FUNCTION cpu_min_nproma 
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  ! The following three functions are for conversion of 1D to 2D indices and vice versa
  !
  ! Treatment of 0 (important for empty patches) and negative numbers:
  !
  ! Converting 1D => 2D:
  !
  ! 0 always is mapped to blk_no = 1, idx_no = 0
  ! negative numbers: Convert usings ABS(j) and negate idx_no
  !
  ! Thus: blk_no >= 1 always!
  !       idx_no > 0  for j > 0
  !       idx_no = 0  for j = 0
  !       idx_no < 0  for j < 0
  !
  ! This mimics mostly the behaviour of reshape_idx in mo_model_domimp_patches
  ! with a difference for nproma=1 and j=0 (where reshape_idx returns blk_no=0, idx_no=1)
  !
  ! The consisten treatment of 0 in the above way is very important for empty patches
  ! where start_index=1, end_index=0
  !
  ! Converting 2D => 1D:
  ! Trying to invert the above and catching cases with blk_no < 1
  !-------------------------------------------------------------------------
  ELEMENTAL INTEGER FUNCTION blk_no(j)
#if defined(_OPENACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(IN) :: j
    blk_no = MAX((ABS(j)-1)/nproma + 1, 1) ! i.e. also 1 for j=0, nproma=1
  END FUNCTION blk_no

  ELEMENTAL INTEGER FUNCTION idx_no(j)
#if defined(_OPENACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(IN) :: j
    IF(j==0) THEN
      idx_no = 0
    ELSE
      idx_no = SIGN(MOD(ABS(j)-1,nproma)+1, j)
    ENDIF
  END FUNCTION idx_no

  ELEMENTAL INTEGER FUNCTION idx_1d(jl,jb)
#if defined(_OPENACC)
    !$ACC ROUTINE SEQ
#endif
    INTEGER, INTENT(IN) :: jl, jb
    IF(jb<=0) THEN
      idx_1d = 0 ! This covers the special case nproma==1,jb=0,jl=1
                 ! All other cases are invalid and get also a 0 returned
    ELSE
      idx_1d = SIGN((jb-1)*nproma + ABS(jl), jl)
    ENDIF
  END FUNCTION idx_1d

END MODULE mo_parallel_config
