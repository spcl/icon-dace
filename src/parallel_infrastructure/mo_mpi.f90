! Basic module initializing the MPI communication and handling most
! of the MPI calls (wrapper functions).
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
!
!  MPI-Handshake
!  -------------
!
!  To allow the co-existence of external processes ICON proceeds the
!  MPI-Handshake (see https://gitlab.dkrz.de/dkrz-sw/mpi-handshake)
!  with the group name "icon" directly after initializing MPI.  This
!  yields a MPI_Comm where all processes are contained in that also
!  have provided this group name.

!  Furthermore ICON participates in the groups "yac" and/or "comin" if
!  these modules are enabled. The resulting communicators are then
!  used to initialize the respective software component.
!
!  In the following the further split-up of the icon communicator is described.
!
!  ICON communicator split
!  -----------------------
!
!  ICON processors are divided into
!    1.    worker PEs    : majority of MPI tasks, doing the actual work
!    2.    I/O PEs       : dedicated I/O server tasks          (only for parallel_nml::num_io_procs > 0)
!    3.    one test PE   : for verification runs               (only for parallel_nml::p_test_run == .TRUE.)
!    4.    restart PEs   : for asynchronous restart writing    (only for dedicatedRestartProcs > 0)
!    5.    prefetch PEs  : for prefetching of data             (only for parallel_nml::num_prefetch_proc > 0)
!    6.    radar I/O PEs : for some async. tasks of radar      (only for parallel_nml::num_io_procs_radar > 0)
!                          forward operator EMVORADO (asynchr. processing: I/O and postprocessing tasks like
!                          radar composites, bubble generator, superobservations)
!
!  The communicators are split like this:
!
!          0    p_work_pe0    p_io_pe0    p_restart_pe0    p_pref_pe0    p_radario_pe0   process_mpi_all_size
!
!          |         |            |             |               |                |              |
!          V         V            V             V               V                V              V
!
!          +---------------------------------------------------------------------+---------------+
!          !                      process_mpi_all_comm                                           !
!          +---------+------------+-------------+---------------+----------------+---------------+
!          | test PE | worker PEs |   I/O PEs   |  restart PEs  |  prefetch PEs  | radar I/O PEs |
!          +---------+------------+-------------+---------------+----------------+---------------+
!          |    A    |     B      |     C       |      D        |       E        |               | p_comm_work
!          |    A    |     A      |             |               |                |               | p_comm_work_test
!          |    A    |     B      |     B       |               |                | (see the      | p_comm_work_io
!          |    A    |            |     B       |               |                |    Note^^     | p_comm_io (B is worker PE 0 if num_io_procs == 0)
!          |         |     A      |             |      A        |                |      below)   | p_comm_work_restart
!          |         |     A      |             |               |       A        |               | p_comm_work_pref
!          +---------+------------+-------------+---------------+----------------+---------------+
!
!  Note that there are actually two different p_comm_work_io communicators:
!  One that spans the worker AND I/O PEs as the NAME implies, AND one that IS ONLY defined on the test PE.
!  Similarly, there IS another ghost communicator p_comm_io defined on the test PE.
!  This has the consequence that `my_process_is_io()` IS NOT equivalent to `p_comm_io /= MPI_COMM_NULL`
!
!  Process groups with specific main procs, all of these are called from mo_atmo_model:
!    * I/O (mo_name_list_output): name_list_io_main_proc()
!    * restart (mo_async_restart): restart_main_proc()
!    * prefetch (mo_async_latbc): prefetch_main_proc()
!
!
!  List of MPI communicators:
!  --------------------------
!
!       global_mpi_communicator
!
!         description  : MPI communicator spanning all PEs running  (= MPI_COMM_WORLD).
!         size         : global_mpi_size
!         this PE's ID : my_global_mpi_id (= get_my_global_mpi_id())
!
!
!       process_mpi_all_comm  (= get_my_mpi_all_communicator())
!
!         description  : MPI communicator containing all PEs that are running this
!                        model component. Different from global_mpi_communicator,
!                        if master_nml::no_of_models > 1
!         size         : process_mpi_all_size
!         this PE's ID : my_process_mpi_all_id (= get_my_mpi_all_id())
!
!
!       p_comm_work  (=get_my_mpi_work_communicator())
!
!         description  : MPI communicator for work group. On I/O and test PEs this
!                        defaults to process_mpi_all_comm.
!         size         : num_work_procs (on worker PEs: num_work_procs==p_n_work)
!         this PE's ID : p_pe_work (= get_my_mpi_work_id())
!
!       In the case of the NEC hybrid mode (detached PE0), the sub-communicator p_comm_work_only encompasses
!       all true worker PEs (excluding PE0), otherwise, p_comm_work_only is identical to p_comm_work
!
!
!       ... less important MPI communicators ...
!
!       p_comm_work_test (size = 0/1)
!         description  : MPI communicator spanning work group and test PE
!                        in verification mode parallel_nml::p_test_run == .TRUE.
!       p_comm_work_io
!         description  : MPI communicator spanning work group and I/O PEs
!       p_comm_io
!         description  : MPI communicator spanning I/O PEs
!       p_comm_work_2_io
!         description  : Inter(!)communicator work PEs - I/O PEs
!       p_comm_work_pref
!         description  : MPI Communicator spanning work group and prefetch PEs
!       p_comm_work_2_pref
!         description  : Inter(!)communicator work PEs - prefetching PEs
!
!
!  Processor splitting:
!  --------------------
!
!       In order to improve the parallel load balancing, the domains
!       of the first refinement level can be distributed to disjoint
!       processor subsets. This distribution process is weighted by
!       the namelist parameter grid_nml::patch_weight.  Since the
!       processor sets are disjoint, MPI calls happen independently on
!       each PE subset.
!
!       For global operations (sum, min, max in mo_sync), each patch
!       contains additional data members which are initialized in
!       mo_complete_subdivision::set_patch_communicators:
!
!       p_patch % comm
!         description  : MPI communicator for this patch. If no processor
!                        splitting enabled: p_patch%comm==p_comm_work.
!         size         : p_patch%n_proc
!         this PE's ID : p_patch%rank
!
!       p_patch % proc0
!         description  : global id of processor with rank 0 (within the
!                        working set p_comm_work)
!
!       The global communicator which is currently in use is stored by
!       a call to mo_mpi::push_glob_comm in the data structure
!       mo_mpi::glob_comm(0:max_lev), where the level number is equal
!       to 0 for the global grid, 1 for the first generation of child
!       domains etc.
!
!
! Split of process_mpi_all_comm and stdio process:
! ------------------------------------------------
!
!    The process_mpi_all_comm is the whole model communicator and is split
!    to test/work/io/restart, in this order
!
!    The process_mpi_stdio_id is always the 0 process of the process_mpi_all_comm
!    If is p_test, then the testing process is also the stdio process
!
!    The my_process_is_mpi_workroot() is true for the 0 process of the
!    p_comm_work communicator, this communicator does not include the test process,
!    the io and the restart processes. This communicator and and the mpi_workroot
!    process should be used in the case of gather and scatter procedures that do
!    not use the io process (the third part). This will be different from using
!    the stdio process only in the case of p_test
!
!
! ^^A Note on radar I/O PEs:
! --------------------------
!
!   If ANY(run_nml::luse_radarfwo(1:max_dom)) and parallel_nml::num_io_procs_radar > 0,
!   These PEs are appended at the very end of the PE list, but are not used by ICON itself.
!   They are provided exclusively to the radar forward operator EMVORADO for
!   Input/Output of observed and simulated radar data and postprocessing tasks
!   such as computing superobservations for data assimilation,
!   generating radar reflectivity composites for model verification and visualization,
!   or the "warm bubble generator" for triggering missing convective cells.
!   From this set of PEs, EMVORADO creates different communicators internally on its own,
!   triggered by a "CALL init_emvorado_mpi()". These communicators
!   are only used in EMVORADO and its interface (src/data_assimilation/interfaces/) and consist of:
!
!     - a clone of p_comm_work                    : icomm_cart_fwo
!     - a communicator for the radar I/O Pes      : icomm_radario
!     - a combined communicator                   : icomm_radar = icomm_cart_fwo + icomm_radario
!
!   In case EMVORADO runs asynchronously for more than 1 (N) different model domains, i.e.,
!   run_nml::luse_radarfwo(1:N) = .TRUE. and parallel_nml::num_io_procs_radar > 0,
!   and if parallel_nml::num_io_procs_radar > N,
!   there are also separate radar I/O sub-communicators for each model domain:
!
!     - a subset of icomm_radario for each domain : icomm_radario_dom(i), i=1...N
!     - a combined communicator for each domain   : icomm_radar_dom(i) = icomm_cart_fwo + icomm_radario_dom(i), i=1...N
!     - icomm_radario is still available as the superset of all icomm_radario_dom(i)
!     - icomm_radar is still available as icomm_cart_fwo + icomm_radario
!
!   This enables to run the radar I/O and postprocessing for each domain
!   separately and in parallel, minimizing communication latency on the worker
!   procs.
!
!   Although ICON does not need to access these communicators explicitly,
!   it checks at a few places if the PE does belong to icomm_radar or icomm_radario.
!   For this purpose, the following is provided PUBLIC below:
!
!   process_mpi_radario_size        : equals parallel_nml::num_io_procs_radar
!   process_mpi_all_radarioroot_id  : ID of root of icomm_radario in the global communicator (= p_radario_pe0)
!   my_process_is_radar()           : returns .TRUE. if PE belongs to icomm_radar
!   my_process_is_radario()         : returns .TRUE. if PE belongs to icomm_radario
!   my_process_is_mpi_radarioroot() : returns .TRUE. if PE is root of icomm_radario
!

!This IS a small helper to avoid a full #ifdef ... #ELSE ... #endif sequence where we
!can just replace the MPI symbol with something constant.






MODULE mo_mpi

  ! Comment: Please use basic WRITE to nerr for messaging in the whole
  !          MPI package to achieve proper output.

  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_signed_char, c_int
  








  USE mo_kind, ONLY: i4, i8, dp, sp, wp
  USE mo_io_units,       ONLY: nerr
  USE mo_impl_constants, ONLY: pio_type_async, pio_type_cdipio
  USE mtime,             ONLY: datetime, max_datetime_str_len, datetimeToString, &
    &                          newDatetime, deallocateDatetime









  USE mo_util_system, ONLY: util_exit

  USE mo_master_control, ONLY: get_my_process_type, hamocc_process, ocean_process, process_exists, &
       &                       my_process_is_hamocc, my_process_is_ocean

  USE mo_coupling, ONLY: init_coupler, finalize_coupler
  USE mo_exception,           ONLY: init_logger


  IMPLICIT NONE

  PRIVATE                          ! all declarations are private


  ! start/stop methods
  PUBLIC :: start_mpi, stop_mpi, abort_mpi


  ! split the global communicator to _process_mpi_communicator
  PUBLIC :: split_global_mpi_communicator
  !The given communicator will be the all communicator for this component
!   PUBLIC :: set_process_mpi_communicator
  ! Sets the test, work, i/o and prefetch communicators
  PUBLIC :: set_mpi_work_communicators
  ! set other parameters
  PUBLIC :: set_process_mpi_name
  ! generates a intercomm between the work PEs of two model components
  PUBLIC :: get_mpi_work_intercomm

  PUBLIC :: push_glob_comm, pop_glob_comm, get_glob_proc0

  ! Logical functions
  PUBLIC :: run_is_global_mpi_parallel
  PUBLIC :: my_process_is_stdio, my_process_is_mpi_parallel, my_process_is_mpi_all_parallel
  PUBLIC :: my_process_is_mpi_seq, my_process_is_mpi_test, my_process_is_mpi_workroot
  PUBLIC :: my_process_is_mpi_ioroot
  PUBLIC :: my_process_is_mpi_all_seq, my_process_is_io
  PUBLIC :: my_process_is_global_root
  PUBLIC :: my_process_is_mpi_prefroot, my_process_is_pref
  PUBLIC :: my_process_is_restart
  PUBLIC :: my_process_is_work, my_process_is_work_only

  ! get parameters
  PUBLIC :: get_my_global_mpi_communicator   ! essentially a copy of MPI_COMM_WORLD
  PUBLIC :: get_my_mpi_all_communicator   ! the communicator for the specific component, ie the process_mpi_all_comm
  PUBLIC :: get_my_mpi_all_comm_size   ! this is the the size of the communicator for the specific component
  PUBLIC :: get_my_mpi_work_communicator   ! the communicator for the workers of this component
  PUBLIC :: get_my_mpi_work_comm_size   ! this is the the size of the workers
  PUBLIC :: get_hamocc_ocean_mpi_communicator 

  PUBLIC :: get_mpi_all_workroot_id, get_my_global_mpi_id, get_my_mpi_all_id
  PUBLIC :: get_mpi_all_ioroot_id
  PUBLIC :: get_my_mpi_work_id, get_mpi_prefroot_id
  PUBLIC :: default_comm_type, null_comm_type


  ! some public communicators
  PUBLIC :: process_mpi_all_comm
  PUBLIC :: process_mpi_all_test_id, process_mpi_all_workroot_id, &
    &       process_mpi_all_ioroot_id, &
    &       process_mpi_all_prefroot_id, p_comm_work_pref_compute_pe0

  PUBLIC :: p_comm_work, p_comm_work_test, p_comm_work_2_test, p_comm_work_only
  PUBLIC :: p_comm_work_2_io, p_comm_work_io, p_comm_dio_io, &
    &       p_comm_io, p_comm_work_pref, p_comm_work_2_pref
  !restart communicators
  PUBLIC :: p_comm_work_2_restart, p_comm_work_restart
  PUBLIC :: p_communicator_a, p_communicator_b, p_communicator_d


  PUBLIC :: process_mpi_io_size, process_mpi_restart_size, process_mpi_pref_size

  PUBLIC :: process_mpi_all_size
  PUBLIC :: process_mpi_radario_size, process_mpi_all_radarioroot_id
  PUBLIC :: my_process_is_radar, my_process_is_radario, my_process_is_mpi_radarioroot

  PUBLIC :: process_mpi_stdio_id
  PUBLIC :: process_mpi_root_id
  PUBLIC :: process_work_io0
  PUBLIC :: process_work_pref0

  ! Main communication methods
  PUBLIC :: global_mpi_barrier
  PUBLIC :: work_mpi_barrier

  PUBLIC :: p_comm_size
  PUBLIC :: p_comm_rank
  PUBLIC :: p_send, p_recv, p_sendrecv, p_bcast, p_barrier
! PUBLIC :: p_bcast_achar
  PUBLIC :: p_get_bcast_role
  PUBLIC :: p_isend, p_irecv, p_wait, p_wait_any,         &
    &       p_irecv_packed, p_send_packed, p_recv_packed, &
    &       p_bcast_packed,                               &
    &       p_pack_int, p_pack_bool, p_pack_real,         &
    &       p_pack_int_1d, p_pack_real_1d,                &
    &       p_pack_string, p_pack_real_2d,                &
    &       p_pack_size_int, p_pack_size_bool,            &
    &       p_pack_size_real_dp, p_pack_size_string,         &
    &       p_unpack_int, p_unpack_bool, p_unpack_real,   &
    &       p_unpack_int_1d, p_unpack_real_1d,            &
    &       p_unpack_string, p_unpack_real_2d, p_test
  PUBLIC :: p_max, p_min, p_lor, p_sum, p_global_sum, p_field_sum
  PUBLIC :: p_scatter, p_gather
  PUBLIC :: p_probe
  PUBLIC :: p_gatherv, p_allgather, p_allgatherv
  PUBLIC :: p_scatterv
  PUBLIC :: p_allreduce_max
  PUBLIC :: p_reduce
  PUBLIC :: p_allreduce
  PUBLIC :: p_commit_type_struct
  PUBLIC :: p_alltoall
  PUBLIC :: p_alltoallv
  PUBLIC :: p_alltoallv_p2p
  PUBLIC :: p_clear_request
  PUBLIC :: p_mpi_wtime
  PUBLIC :: p_comm_is_intercomm
  PUBLIC :: p_comm_remote_size
  PUBLIC :: get_mpi_comm_world_ranks

  PUBLIC :: p_isEqual

  !----------- to be removed -----------------------------------------
  PUBLIC :: p_pe, p_io
  PUBLIC :: num_test_procs, num_work_procs,     &
       &    p_work_pe0, p_io_pe0,               &
       &    p_n_work, p_pe_work, &
       &    p_pref_pe0
  PUBLIC :: p_pe_work_only
  !--------------------------------------------------------------------

  PUBLIC :: MPI_ANY_SOURCE, MPI_COMM_NULL, MPI_COMM_SELF
  PUBLIC :: MPI_REQUEST_NULL

  ! real data type matching real type of MPI implementation
  PUBLIC :: p_real_dp, p_real_sp, p_real
  PUBLIC :: p_int
  PUBLIC :: p_int_i4
  PUBLIC :: p_int_i8
  PUBLIC :: p_bool
  PUBLIC :: p_address_kind
  PUBLIC :: p_real_dp_byte, p_real_sp_byte, p_int_byte
  PUBLIC :: p_int_i4_byte, p_int_i8_byte
  PUBLIC :: p_mpi_comm_null
  
  ! mpi reduction operators
  PUBLIC :: mpi_lor, mpi_land, mpi_sum, mpi_min, mpi_max, &
       mpi_minloc, mpi_maxloc
  INTEGER, PARAMETER :: mpi_lor = 1, mpi_land = 2, mpi_sum = 3, &
       mpi_min = 4, mpi_max = 5, mpi_minloc = 6, mpi_maxloc = 7

  ! NCCL-related routines. Typically no-op when NCCL is not active.
  PUBLIC :: push_gpu_comm_queue, pop_gpu_comm_queue, get_comm_acc_queue, &
            acc_wait_comms, comm_group_start, comm_group_end

  PUBLIC ::get_mpi_time,ellapsed_mpi_time

  TYPE t_work_root_process
    INTEGER :: comp_id
    INTEGER :: global_mpi_id
  END TYPE t_work_root_process

  TYPE (t_work_root_process), ALLOCATABLE :: p_work_root_processes(:)

  ! general run time information

  INTEGER, PARAMETER :: p_address_kind = i8    ! should not get touched at all
  INTEGER, PARAMETER :: MPI_COMM_NULL  = 0
  INTEGER, PARAMETER :: MPI_COMM_SELF  = 1
  ! dummy arguments for function calls:
  INTEGER, PARAMETER :: MPI_ANY_SOURCE = 0
  INTEGER, PARAMETER :: mpi_request_null = 0
  ! this is the global communicator
  INTEGER, PARAMETER :: global_mpi_communicator = 0  ! replaces MPI_COMM_WORLD

  INTEGER :: hamocc_ocean_mpi_communicator

  ! public parallel run information
  CHARACTER(len=64) :: global_mpi_name
  CHARACTER(len=64) :: process_mpi_name
  ! Do not change the following parameters
  INTEGER, PARAMETER :: process_mpi_stdio_id = 0
  INTEGER, PARAMETER :: process_mpi_root_id = 0
  INTEGER, PARAMETER :: default_comm_type = 1
  INTEGER, PARAMETER :: null_comm_type = 0

  ! communicator sets
  INTEGER :: global_mpi_size          ! total number of processes in global world
  INTEGER :: my_global_mpi_id         ! process id in global world
  LOGICAL :: is_global_mpi_parallel

  ! this is the communicator for the whole component (atmo/ocean/etc)
  INTEGER :: process_mpi_all_comm     ! communicator in the whole model-component
  INTEGER :: process_mpi_all_size     ! total number of processes in the whole model-component
  INTEGER :: my_process_mpi_all_id        ! the id in the whole model communicator
  INTEGER :: process_mpi_all_workroot_id  ! the root process in component
  INTEGER :: process_mpi_all_ioroot_id    ! the first I/O process
  INTEGER :: process_mpi_all_radarioroot_id    ! the first radar I/O process
  INTEGER :: process_mpi_all_test_id  ! the test process in component
  INTEGER :: process_work_io0
  INTEGER :: process_mpi_all_prefroot_id ! the id of the first prefetch process
  INTEGER :: p_comm_work_pref_compute_pe0 ! the ID for Communicator spanning work group and prefetch PEs
  INTEGER :: process_work_pref0           ! ID of prefetch PE0 within p_comm_work_pref communicator

  LOGICAL :: process_is_mpi_parallel
  LOGICAL :: process_is_stdio

  ! this is the local work communicator (computation, i/o, etc)
!   INTEGER :: process_mpi_local_comm     ! communicator in the work group
!   INTEGER :: process_mpi_local_size     ! total number of processes in the whole model-component
!   INTEGER :: my_process_mpi_local_id

  INTEGER :: my_mpi_function  ! test, work, i/o, restart_output or prefetch
  INTEGER, PARAMETER :: test_mpi_process = 1
  INTEGER, PARAMETER :: work_mpi_process = 2
  INTEGER, PARAMETER :: io_mpi_process = 3
  INTEGER, PARAMETER :: restart_mpi_process = 4
  INTEGER, PARAMETER :: pref_mpi_process = 5
  INTEGER, PARAMETER :: radario_mpi_process = 6

  !------------------------------------------------------------
  ! Processor distribution:
  ! num_test_procs:      0 or 1
  ! num_work_procs:      number of procs running in parallel on the model
  ! process_mpi_io_size: number of procs for I/O
  ! num_restart_procs:   number of procs used for writing restart files
  ! num_prefetch_proc:      number of procs used for prefetching of data
  ! num_test_procs + num_work_procs + process_mpi_io_size + num_restart_procs + num_prefetch_proc = process_mpi_all_size
  INTEGER :: num_test_procs
  INTEGER :: num_work_procs
  INTEGER :: p_radario_pe0
  INTEGER :: process_mpi_radario_size
  INTEGER :: process_mpi_io_size
  INTEGER :: process_mpi_restart_size
  INTEGER :: process_mpi_pref_size

  ! Note: p_work_pe0, p_io_pe0 are identical on all PEs

  INTEGER :: p_work_pe0    ! Number of workgroup PE 0 within all PEs
  INTEGER :: p_workonly_pe0 ! Number of workgroup PE 0 excluding detached PE0
  INTEGER :: p_io_pe0      ! Number of I/O PE 0 within all PEs (process_mpi_all_size if no I/O PEs)
  INTEGER :: p_pref_pe0    ! Number of the prefetching PE 0 within all PEs (process_mpi_all_size
  !                          if no prefetching PEs)

  ! Note: p_n_work, p_pe_work are NOT identical on all PEs

  ! p_n_work: Number of PEs working together:
  ! - num_work_procs    for non-verification runs
  ! - num_work_procs    for verification runs on pes != p_test_pe
  ! - num_test_procs    for verification runs on p_test_pe
  ! - num_io_procs      always on I/O pes
  ! - num_restart_procs on restart PEs

  INTEGER :: p_n_work
  INTEGER :: p_pe_work        ! PE number within work group
  INTEGER :: p_pe_work_only   ! PE number within work group excluding detached stdio-PE (NEC hybrid mode), otherwise same as p_pe_work


  ! MPI communicators
  INTEGER :: p_comm_work           ! Communicator for work group
  INTEGER :: p_comm_work_only      ! Communicator for work group excluding detached stdio-PE (NEC hybrid mode), otherwise same as p_comm_work
  INTEGER :: p_comm_work_test      ! Communicator spanning work group and test PE
  INTEGER :: p_comm_work_2_test
  INTEGER :: p_comm_work_io        ! Communicator spanning work group and I/O PEs
  INTEGER :: p_comm_dio_io         ! Communicator spanning detached std I/O and I/O PEs
  INTEGER :: p_comm_io             ! Communicator spanning the I/O PEs
  INTEGER :: p_comm_work_2_io      ! Inter(!)communicator work PEs - I/O PEs
  INTEGER :: p_comm_work_restart   ! Communicator spanning work group and Restart Output PEs
  INTEGER :: p_comm_work_2_restart ! Inter(!)communicator work PEs - Restart PEs
  INTEGER :: p_comm_work_pref           ! Communicator spanning work group and prefetch PEs
  INTEGER :: p_comm_work_2_pref    ! Inter(!)communicator work PEs - prefetching PEs

  INTEGER :: p_communicator_a ! for Set A
  INTEGER :: p_communicator_b ! for Set B
  INTEGER :: p_communicator_d ! for debug node

  INTEGER :: p_pe     = 0     ! this is the PE number of this task
  INTEGER :: p_io     = 0     ! PE number of PE handling IO


! non blocking calls

  ! module intrinsic names

!!#ifndef 1
!!LK  INTEGER :: iope                  ! PE able to do IO
!!#endif


  ! MPI native transfer types

  INTEGER :: p_int     = 0
  INTEGER :: p_real    = 0
  INTEGER :: p_bool    = 0
  INTEGER :: p_char    = 0

  ! MPI transfer types

  INTEGER :: p_real_dp = 0
  INTEGER :: p_real_sp = 0
  INTEGER :: p_int_i4  = 0
  INTEGER :: p_int_i8  = 0

  ! MPI native transfer types byte size

  INTEGER :: p_int_byte     = 0
  INTEGER :: p_real_byte    = 0
!  INTEGER :: p_bool_byte    = 0
!  INTEGER :: p_char_byte    = 0

  ! MPI transfer types byte size

  INTEGER :: p_real_dp_byte = 0
  INTEGER :: p_real_sp_byte = 0
  INTEGER :: p_int_i4_byte  = 0
  INTEGER :: p_int_i8_byte  = 0

  INTEGER :: p_mpi_comm_null = -32766
  
  ! Flag if processor splitting is active
  LOGICAL, PUBLIC :: proc_split = .FALSE.
  LOGICAL, PUBLIC :: i_am_accel_node = .FALSE.

  LOGICAL, PARAMETER :: nccl_active = .FALSE.


  ! communicator stack for global sums
  INTEGER, PARAMETER :: max_lev = 10 ! 2 is sufficient
  INTEGER, PUBLIC :: comm_lev = 0, glob_comm(0:max_lev), comm_proc0(0:max_lev)

  ! Storage for buffered non-blocking point-to-point communication. This is
  ! especially needed for communication between I/O servers. These tasks may
  ! be completely asynchronous and therefore ISENDs may be launched before a
  ! corresponding IRECVs has been issued.
  INTEGER :: mpi_buffer(10000)


  ! define generic interfaces to allow proper compiling with picky compilers
  ! like NAG f95 for clean argument checking and shortening the call sequence.

  INTERFACE p_send
     MODULE PROCEDURE p_send_char
     MODULE PROCEDURE p_send_real
     MODULE PROCEDURE p_send_sreal
     MODULE PROCEDURE p_send_int
     MODULE PROCEDURE p_send_bool
     MODULE PROCEDURE p_send_real_1d
     MODULE PROCEDURE p_send_sreal_1d
     MODULE PROCEDURE p_send_int_1d
     MODULE PROCEDURE p_send_bool_1d
     MODULE PROCEDURE p_send_char_1d
     MODULE PROCEDURE p_send_real_2d
     MODULE PROCEDURE p_send_int_2d
     MODULE PROCEDURE p_send_bool_2d
     MODULE PROCEDURE p_send_real_3d
     MODULE PROCEDURE p_send_int_3d
     MODULE PROCEDURE p_send_bool_3d
     MODULE PROCEDURE p_send_real_4d
     MODULE PROCEDURE p_send_int_4d
     MODULE PROCEDURE p_send_bool_4d
     MODULE PROCEDURE p_send_real_5d
  END INTERFACE

  INTERFACE p_isend
     MODULE PROCEDURE p_isend_char
     MODULE PROCEDURE p_isend_real
     MODULE PROCEDURE p_isend_sreal
     MODULE PROCEDURE p_isend_int
     MODULE PROCEDURE p_isend_bool
     MODULE PROCEDURE p_isend_char_1d
     MODULE PROCEDURE p_isend_real_1d
     MODULE PROCEDURE p_isend_sreal_1d
     MODULE PROCEDURE p_isend_int_1d
     MODULE PROCEDURE p_isend_bool_1d
     MODULE PROCEDURE p_isend_real_2d
     MODULE PROCEDURE p_isend_sreal_2d
     MODULE PROCEDURE p_isend_int_2d
     MODULE PROCEDURE p_isend_bool_2d
     MODULE PROCEDURE p_isend_real_3d
     MODULE PROCEDURE p_isend_int_3d
     MODULE PROCEDURE p_isend_bool_3d
     MODULE PROCEDURE p_isend_real_4d
     MODULE PROCEDURE p_isend_int_4d
     MODULE PROCEDURE p_isend_bool_4d
     MODULE PROCEDURE p_isend_real_5d
  END INTERFACE

  INTERFACE p_recv
     MODULE PROCEDURE p_recv_char
     MODULE PROCEDURE p_recv_real
     MODULE PROCEDURE p_recv_sreal
     MODULE PROCEDURE p_recv_int
     MODULE PROCEDURE p_recv_bool
     MODULE PROCEDURE p_recv_real_1d
     MODULE PROCEDURE p_recv_sreal_1d
     MODULE PROCEDURE p_recv_int_1d
     MODULE PROCEDURE p_recv_bool_1d
     MODULE PROCEDURE p_recv_char_1d
     MODULE PROCEDURE p_recv_real_2d
     MODULE PROCEDURE p_recv_int_2d
     MODULE PROCEDURE p_recv_bool_2d
     MODULE PROCEDURE p_recv_real_3d
     MODULE PROCEDURE p_recv_int_3d
     MODULE PROCEDURE p_recv_bool_3d
     MODULE PROCEDURE p_recv_real_4d
     MODULE PROCEDURE p_recv_int_4d
     MODULE PROCEDURE p_recv_bool_4d
     MODULE PROCEDURE p_recv_real_5d
  END INTERFACE

  INTERFACE p_irecv
     MODULE PROCEDURE p_irecv_char
     MODULE PROCEDURE p_irecv_real
     MODULE PROCEDURE p_irecv_sreal
     MODULE PROCEDURE p_irecv_int
     MODULE PROCEDURE p_irecv_bool
     MODULE PROCEDURE p_irecv_char_1d
     MODULE PROCEDURE p_irecv_real_1d
     MODULE PROCEDURE p_irecv_sreal_1d
     MODULE PROCEDURE p_irecv_int_1d
     MODULE PROCEDURE p_irecv_bool_1d
     MODULE PROCEDURE p_irecv_real_2d
     MODULE PROCEDURE p_irecv_sreal_2d
     MODULE PROCEDURE p_irecv_int_2d
     MODULE PROCEDURE p_irecv_bool_2d
     MODULE PROCEDURE p_irecv_real_3d
     MODULE PROCEDURE p_irecv_int_3d
     MODULE PROCEDURE p_irecv_bool_3d
     MODULE PROCEDURE p_irecv_real_4d
     MODULE PROCEDURE p_irecv_int_4d
     MODULE PROCEDURE p_irecv_bool_4d
  END INTERFACE p_irecv

  INTERFACE p_wait
    MODULE PROCEDURE p_wait
    MODULE PROCEDURE p_wait_1
    MODULE PROCEDURE p_wait_n
    MODULE PROCEDURE p_wait_n_stati
  END INTERFACE p_wait

  INTERFACE p_clear_request
    MODULE PROCEDURE p_clear_request
    MODULE PROCEDURE p_clear_requests
  END INTERFACE p_clear_request

  INTERFACE p_sendrecv
     MODULE PROCEDURE p_sendrecv_real_1d
     MODULE PROCEDURE p_sendrecv_real_2d
     MODULE PROCEDURE p_sendrecv_real_3d
     MODULE PROCEDURE p_sendrecv_real_4d
     MODULE PROCEDURE p_sendrecv_char_array
  END INTERFACE

  INTERFACE p_send_packed
    MODULE PROCEDURE p_send_packed
    MODULE PROCEDURE p_send_packed_2d
  END INTERFACE p_send_packed

  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_real
     MODULE PROCEDURE p_bcast_real_single
     MODULE PROCEDURE p_bcast_int_i4
     MODULE PROCEDURE p_bcast_int_i8
     MODULE PROCEDURE p_bcast_bool
     MODULE PROCEDURE p_bcast_real_1d
     MODULE PROCEDURE p_bcast_real_1d_single
     MODULE PROCEDURE p_bcast_int_1d
     MODULE PROCEDURE p_bcast_int_i8_1d
     MODULE PROCEDURE p_bcast_bool_1d
     MODULE PROCEDURE p_bcast_real_2d
     MODULE PROCEDURE p_bcast_real_2d_single
     MODULE PROCEDURE p_bcast_int_2d
     MODULE PROCEDURE p_bcast_bool_2d
     MODULE PROCEDURE p_bcast_real_3d
     MODULE PROCEDURE p_bcast_int_3d
     MODULE PROCEDURE p_bcast_bool_3d
     MODULE PROCEDURE p_bcast_real_4d
     MODULE PROCEDURE p_bcast_int_4d
     MODULE PROCEDURE p_bcast_bool_4d
     MODULE PROCEDURE p_bcast_real_5d
     MODULE PROCEDURE p_bcast_int_7d
     MODULE PROCEDURE p_bcast_real_7d
     MODULE PROCEDURE p_bcast_char
     MODULE PROCEDURE p_bcast_cchar
     MODULE PROCEDURE p_bcast_char_1d
     MODULE PROCEDURE p_bcast_datetime
     MODULE PROCEDURE p_bcast_char_2d
  END INTERFACE

  INTERFACE p_scatter
     MODULE PROCEDURE p_scatter_real_1d1d
     MODULE PROCEDURE p_scatter_real_2d1d
     MODULE PROCEDURE p_scatter_sp_1d1d
     MODULE PROCEDURE p_scatter_sp_2d1d
     MODULE PROCEDURE p_scatter_int_1d1d
     MODULE PROCEDURE p_scatter_int_2d1d
  END INTERFACE

  INTERFACE p_gather
     MODULE PROCEDURE p_gather_real_0d1d
     MODULE PROCEDURE p_gather_real_1d2d
     MODULE PROCEDURE p_gather_real_2d3d
     MODULE PROCEDURE p_gather_real_5d6d
     MODULE PROCEDURE p_gather_real_1d1d
     MODULE PROCEDURE p_gather_int_0d1d
     MODULE PROCEDURE p_gather_int_1d1d
     MODULE PROCEDURE p_gather_int_1d2d
     MODULE PROCEDURE p_gather_int_2d3d
     MODULE PROCEDURE p_gather_char_0d1d
     MODULE PROCEDURE p_gather_bool_0d1d
  END INTERFACE

  INTERFACE p_allgather
     MODULE PROCEDURE p_allgather_int_0d1d
     MODULE PROCEDURE p_allgather_int_1d2d
  END INTERFACE

  INTERFACE p_allgatherv
     MODULE PROCEDURE p_allgatherv_real_1d
     MODULE PROCEDURE p_allgatherv_int_1d
     MODULE PROCEDURE p_allgatherv_int_1d_contiguous
  END INTERFACE

  INTERFACE p_scatterv
    MODULE PROCEDURE p_scatterv_real1D2D
    MODULE PROCEDURE p_scatterv_real1D1D
    MODULE PROCEDURE p_scatterv_int_1d1d
    MODULE PROCEDURE p_scatterv_single1D1D
  END INTERFACE

  INTERFACE p_gatherv
    MODULE PROCEDURE p_gatherv_int
    MODULE PROCEDURE p_gatherv_real2D1D
    MODULE PROCEDURE p_gatherv_real3D1D
    MODULE PROCEDURE p_gatherv_int2D1D
    MODULE PROCEDURE p_gatherv_real2D2D
    MODULE PROCEDURE p_gatherv_sreal2D2D
    MODULE PROCEDURE p_gatherv_int2D2D
  END INTERFACE

  INTERFACE p_max
     MODULE PROCEDURE p_max_0d
     MODULE PROCEDURE p_max_int_0d
     MODULE PROCEDURE p_max_1d
     MODULE PROCEDURE p_max_int_1d
     MODULE PROCEDURE p_max_2d
     MODULE PROCEDURE p_max_3d
     MODULE PROCEDURE p_max_0d_sp
     MODULE PROCEDURE p_max_1d_sp
     MODULE PROCEDURE p_max_2d_sp
     MODULE PROCEDURE p_max_3d_sp
  END INTERFACE

  INTERFACE p_min
     MODULE PROCEDURE p_min_0d
     MODULE PROCEDURE p_min_int_0d
     MODULE PROCEDURE p_min_1d
     MODULE PROCEDURE p_min_int_1d
     MODULE PROCEDURE p_min_2d
     MODULE PROCEDURE p_min_3d
  END INTERFACE

  INTERFACE p_lor
    MODULE PROCEDURE p_lor_0d
  END INTERFACE p_lor

  INTERFACE p_sum
     MODULE PROCEDURE p_sum_sp_0d
     MODULE PROCEDURE p_sum_sp_1d
     MODULE PROCEDURE p_sum_dp_0d
     MODULE PROCEDURE p_sum_dp_1d
     MODULE PROCEDURE p_sum_dp_2d
     MODULE PROCEDURE p_sum_dp_3d
     MODULE PROCEDURE p_sum_i8_1d
     MODULE PROCEDURE p_sum_i_1d
     MODULE PROCEDURE p_sum_i_0d
  END INTERFACE

  INTERFACE p_global_sum
     MODULE PROCEDURE p_global_sum_1d
  END INTERFACE

  INTERFACE p_field_sum
     MODULE PROCEDURE p_field_sum_1d
     MODULE PROCEDURE p_field_sum_2d
  END INTERFACE

  !> generic interface for MPI communication calls
  INTERFACE p_allreduce_max
    MODULE PROCEDURE p_allreduce_max_int_1d
  END INTERFACE

  INTERFACE p_reduce
    MODULE PROCEDURE p_reduce_i8_0d
    MODULE PROCEDURE p_reduce_i4_0d
  END INTERFACE p_reduce

  INTERFACE p_allreduce
    MODULE PROCEDURE p_allreduce_bool_0d
    MODULE PROCEDURE p_allreduce_int4_0d
    MODULE PROCEDURE p_allreduce_int8_0d
  END INTERFACE p_allreduce

  INTERFACE p_alltoall
    MODULE PROCEDURE p_alltoall_int
  END INTERFACE

  INTERFACE p_alltoallv
    MODULE PROCEDURE p_alltoallv_int
    MODULE PROCEDURE p_alltoallv_int_i8_1d
    MODULE PROCEDURE p_alltoallv_real_2d
    MODULE PROCEDURE p_alltoallv_sreal_2d
    MODULE PROCEDURE p_alltoallv_int_2d
  END INTERFACE

  INTERFACE p_alltoallv_p2p
    MODULE PROCEDURE p_alltoallv_p2p_real_2d
    MODULE PROCEDURE p_alltoallv_p2p_int_2d
  END INTERFACE

  INTERFACE p_isEqual
    MODULE PROCEDURE p_isEqual_int
    MODULE PROCEDURE p_isEqual_charArray
  END INTERFACE

  INTERFACE p_minmax_common
    MODULE PROCEDURE p_minmax_common
    MODULE PROCEDURE p_minmax_common_sp
  END INTERFACE p_minmax_common

  CHARACTER(*), PARAMETER :: modname = 'mo_mpi'

  CHARACTER(len=256) :: message_text = ''
  

CONTAINS


  !------------------------------------------------------------------------------
  ! NCCL routines
  !------------------------------------------------------------------------------


  SUBROUTINE push_gpu_comm_queue(queue)
    INTEGER, INTENT(IN) :: queue
    ! no-op
  END SUBROUTINE

  SUBROUTINE pop_gpu_comm_queue()
    ! no-op
  END SUBROUTINE

  FUNCTION get_comm_acc_queue()
    INTEGER :: get_comm_acc_queue
    get_comm_acc_queue = 1
  END FUNCTION

  SUBROUTINE acc_wait_comms(queue)
    INTEGER, INTENT(IN) :: queue
    !$ACC WAIT(queue)
  END SUBROUTINE

  SUBROUTINE comm_group_start()
    ! no-op
  END SUBROUTINE

  SUBROUTINE comm_group_end()
    ! no-op
  END SUBROUTINE


  !------------------------------------------------------------------------------
  REAL(dp) FUNCTION get_mpi_time()
    get_mpi_time =  0.0_dp
  END FUNCTION get_mpi_time
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  REAL(dp) FUNCTION ellapsed_mpi_time(start_time)
    REAL(dp), INTENT(in) :: start_time
    ellapsed_mpi_time =  0.0_dp
  END FUNCTION ellapsed_mpi_time
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE set_process_mpi_name(name)
    CHARACTER(len=*), INTENT(in) ::name

    process_mpi_name = name

  END SUBROUTINE set_process_mpi_name
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Rene: this has become obsolete and should be erased.
  ! It is used to idendify globally the proccess that caused an error
  INTEGER FUNCTION get_my_global_mpi_id()
    get_my_global_mpi_id = my_global_mpi_id
  END FUNCTION get_my_global_mpi_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_global_mpi_communicator()
    get_my_global_mpi_communicator = global_mpi_communicator
  END FUNCTION get_my_global_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_hamocc_ocean_mpi_communicator()
    get_hamocc_ocean_mpi_communicator = hamocc_ocean_mpi_communicator
  END FUNCTION get_hamocc_ocean_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_communicator()
    get_my_mpi_all_communicator = process_mpi_all_comm
  END FUNCTION get_my_mpi_all_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_comm_size()
    get_my_mpi_all_comm_size = process_mpi_all_size
  END FUNCTION get_my_mpi_all_comm_size
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_all_id()
    get_my_mpi_all_id = my_process_mpi_all_id
  END FUNCTION get_my_mpi_all_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_all_workroot_id()
    get_mpi_all_workroot_id = process_mpi_all_workroot_id
  END FUNCTION get_mpi_all_workroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_all_ioroot_id()
    get_mpi_all_ioroot_id = process_mpi_all_ioroot_id
  END FUNCTION get_mpi_all_ioroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_mpi_prefroot_id()
    get_mpi_prefroot_id = process_mpi_all_prefroot_id
  END FUNCTION get_mpi_prefroot_id
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_communicator()
    get_my_mpi_work_communicator = p_comm_work
  END FUNCTION get_my_mpi_work_communicator
  !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_comm_size()
    get_my_mpi_work_comm_size = p_n_work
  END FUNCTION get_my_mpi_work_comm_size
  !------------------------------------------------------------------------------
  !------------------------------------------------------------------------------
  INTEGER FUNCTION get_my_mpi_work_id()
    get_my_mpi_work_id = p_pe_work
  END FUNCTION get_my_mpi_work_id
  !------------------------------------------------------------------------------



  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_stdio()
    my_process_is_stdio =  .TRUE.
  END FUNCTION my_process_is_stdio
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_io()
    my_process_is_io = (my_mpi_function == io_mpi_process)
  END FUNCTION my_process_is_io
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_radar()
    my_process_is_radar = ( my_process_is_work() .OR. my_process_is_radario() )
  END FUNCTION my_process_is_radar
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_radario()
    my_process_is_radario = (my_mpi_function == radario_mpi_process)
  END FUNCTION my_process_is_radario
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_mpi_radarioroot()
    my_process_is_mpi_radarioroot = (my_process_is_radario() .AND. &
     &                   (my_process_mpi_all_id == process_mpi_all_radarioroot_id))
  END FUNCTION my_process_is_mpi_radarioroot
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_ioroot()
    my_process_is_mpi_ioroot = (my_process_mpi_all_id == process_mpi_all_ioroot_id)
  END FUNCTION my_process_is_mpi_ioroot
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_global_root()
    my_process_is_global_root = (my_global_mpi_id == process_mpi_stdio_id) ! == 0
  END FUNCTION my_process_is_global_root
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_restart()
     my_process_is_restart = (my_mpi_function == restart_mpi_process)
  END FUNCTION my_process_is_restart
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_work()
     my_process_is_work = (my_mpi_function == work_mpi_process)
  END FUNCTION my_process_is_work
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_work_only()
     my_process_is_work_only = (my_process_is_work() .AND. p_pe >= p_workonly_pe0)
  END FUNCTION my_process_is_work_only
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_pref()
    my_process_is_pref = (my_mpi_function == pref_mpi_process)
  END FUNCTION my_process_is_pref
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_prefroot()
    my_process_is_mpi_prefroot = (my_process_is_pref() .AND. &
     &                   (my_process_mpi_all_id == process_mpi_all_prefroot_id))
  END FUNCTION my_process_is_mpi_prefroot
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  LOGICAL FUNCTION my_process_is_mpi_test()
    my_process_is_mpi_test = (my_mpi_function == test_mpi_process)
  END FUNCTION my_process_is_mpi_test
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  ! If is mpi parallel and not a test process
  !! Note: mpi i/o processes do not count in mpi parallel work
  !!        Only computational processes are checked for running in parallel
  LOGICAL FUNCTION my_process_is_mpi_parallel()
    my_process_is_mpi_parallel = process_is_mpi_parallel
  END FUNCTION my_process_is_mpi_parallel
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !>
  ! If is mpi all parallel
  LOGICAL FUNCTION my_process_is_mpi_all_parallel()
    my_process_is_mpi_all_parallel = (process_mpi_all_size > 1)
  END FUNCTION my_process_is_mpi_all_parallel
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  ! If is mpi all sequential
  LOGICAL FUNCTION my_process_is_mpi_all_seq()
    my_process_is_mpi_all_seq = (process_mpi_all_size <= 1)
  END FUNCTION my_process_is_mpi_all_seq
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! If is not mpi work parallel or this is a test process
  !! returns true
  !! Note: mpi i/o processes do not count is mpi parallel work
  !!        Only computational processes are checked for running in parallel
  LOGICAL FUNCTION my_process_is_mpi_seq()
    my_process_is_mpi_seq = .NOT. process_is_mpi_parallel
  END FUNCTION my_process_is_mpi_seq
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  LOGICAL FUNCTION my_process_is_mpi_workroot()
!     my_process_is_mpi_workroot = (my_process_mpi_all_id == process_mpi_all_workroot_id)
    my_process_is_mpi_workroot = (p_pe_work == process_mpi_root_id)
  END FUNCTION my_process_is_mpi_workroot
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  !>
  ! If is not mpi parallel or is a test process
  LOGICAL FUNCTION run_is_global_mpi_parallel()
    run_is_global_mpi_parallel = is_global_mpi_parallel
  END FUNCTION run_is_global_mpi_parallel
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE print_info_stderr (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    IF (my_process_is_stdio() ) THEN
      WRITE (nerr,'(a,a,a,a)') " ", TRIM(name), ": ", TRIM(text)
    ENDIF

  END SUBROUTINE print_info_stderr
  !------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
  SUBROUTINE finish (name, text)
    CHARACTER (len=*), INTENT(in) :: name
    CHARACTER (len=*), INTENT(in) :: text

    WRITE (nerr,'(i4.4,a,a,a,a)') my_global_mpi_id, ": ", TRIM(name), ": ", TRIM(text)
    CALL abort_mpi

  END SUBROUTINE finish
  !------------------------------------------------------------------------------


  !-------------------------------------------------------------------------
  !
  !> Pushes the communicator and proc0 onto the communicator stack.
  !! The communicator stack is needed for global sums if the processor
  !! set is split among different 1st level patches.

  SUBROUTINE push_glob_comm(comm, proc0)
    INTEGER, INTENT(IN) :: comm, proc0

    ! Safety check
    IF(comm_lev>=max_lev) &
      CALL finish('push_glob_comm','max_lev exceeded')

    comm_lev = comm_lev+1
    glob_comm(comm_lev) = comm
    comm_proc0(comm_lev) = proc0

  END SUBROUTINE push_glob_comm

  !-------------------------------------------------------------------------
  !
  !> Pops one level of the communicator stack

  SUBROUTINE pop_glob_comm()

    ! Safety check
    IF(comm_lev<=0) &
      CALL finish('pop_glob_comm','stack empty')

    comm_lev = comm_lev-1

  END SUBROUTINE pop_glob_comm


  !-------------------------------------------------------------------------
  !
  !> @return current top of communicator stack

  FUNCTION get_glob_proc0()
    INTEGER :: get_glob_proc0
    get_glob_proc0 = comm_proc0(comm_lev)
  END FUNCTION get_glob_proc0


  !------------------------------------------------------------------------------
  !-------------------------------------------------------------------------
  ! Warning: The dummy argument num_restart_procs IS NOT identical to the namelist PARAMETER anymore,
  !          rather, it IS the number of *dedicated* restart processes.
  !          We don't care about work processes here that are reused as restart writing processes.
  SUBROUTINE set_mpi_work_communicators(p_test_run, l_test_openmp, num_io_procs,               &
    &                                   num_restart_procs, my_comp_id, num_prefetch_proc,      &
    &                                   num_test_pe, pio_type,                                 &
    &                                   num_io_procs_radar, radar_flag_doms_model, num_dio_procs)

    LOGICAL, INTENT(INOUT)           :: p_test_run, l_test_openmp
    INTEGER, INTENT(INOUT)           :: num_io_procs
    INTEGER, INTENT(INOUT)           :: num_restart_procs
    INTEGER, INTENT(IN)              :: my_comp_id
    INTEGER, INTENT(IN),    OPTIONAL :: num_prefetch_proc, num_test_pe
    INTEGER, INTENT(IN),    OPTIONAL :: pio_type
    INTEGER, INTENT(INOUT), OPTIONAL :: num_io_procs_radar
    LOGICAL, INTENT(IN),    OPTIONAL :: radar_flag_doms_model(:)
    INTEGER, INTENT(IN),    OPTIONAL :: num_dio_procs

!   !local variables
    INTEGER :: my_color, remote_leader, peer_comm, p_error, global_dup_comm
    INTEGER :: tmp_common_intracom, tmp_work_intracom, my_common_intracom_mpi_id, my_work_intracom_mpi_id
    INTEGER :: hamocc_root, ocean_root
    LOGICAL :: I_am_receiver, I_am_sender
    INTEGER :: my_function_comm
    CHARACTER(*), PARAMETER :: method_name = "set_mpi_work_communicators"
    INTEGER :: grp_process_mpi_all_comm, grp_comm_work_io, input_ranks(1), &
               translated_ranks(1), grp_comm_work_pref
    INTEGER :: sizeof_prefetch_processes, pio_type_

    INTEGER :: num_radario_procs
    LOGICAL :: lhave_radar


    INTEGER :: p_restart_pe0 ! Number of Restart Output PE 0 within all PEs
                             ! (process_mpi_all_size if no restart PEs)
    INTEGER :: num_component, i
    INTEGER, ALLOCATABLE :: root_buffer(:)
    INTEGER :: comp_id
    CHARACTER(len=1000) :: message_text = ''
    CHARACTER(len=*), PARAMETER :: &
         routine = modname//'::set_mpi_work_communicators'

    INTEGER :: my_arch, p
    INTEGER, ALLOCATABLE :: glb_arch(:)
    LOGICAL :: arch_mismatch

    IF (PRESENT(num_prefetch_proc)) THEN
      sizeof_prefetch_processes = num_prefetch_proc
    ELSE
      sizeof_prefetch_processes = 0
    ENDIF

    comp_id = my_comp_id

    IF (PRESENT(pio_type)) THEN
      pio_type_ = pio_type
    ELSE
      pio_type_ = pio_type_async
    END IF

    IF (PRESENT(num_io_procs_radar) .AND. PRESENT(radar_flag_doms_model)) THEN
      lhave_radar = .TRUE.
      IF (.NOT.ANY(radar_flag_doms_model)) THEN
        CALL print_info_stderr(method_name, &
             & 'num_io_procs_radar has no effect if there are no radar-active domains!')
        CALL print_info_stderr(method_name, &
             & '--> num_io_procs_radar set to 0')
        num_io_procs_radar = 0
      END IF
    ELSE
      lhave_radar = .FALSE.
    END IF

! check l_test_openmp
    IF (l_test_openmp) THEN
      CALL print_info_stderr(method_name, &
        & 'l_test_openmp has no effect if the model is compiled without OpenMP support')
      CALL print_info_stderr(method_name, &
        & '--> l_test_openmp set to .FALSE.')
      l_test_openmp = .FALSE.
    END IF

    ! check p_test_run, num_io_procs and num_restart_procs
    ! Unconditionally set p_test_run to .FALSE. and num_io_procs to 0,
    ! all other variables are already set correctly
    IF (p_test_run) THEN
      CALL print_info_stderr(method_name, &
       & 'p_test_run has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> p_test_run set to .FALSE.')
      p_test_run = .FALSE.
    END IF
    IF (num_io_procs /= 0) THEN
      CALL print_info_stderr(method_name, &
       & 'num_io_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> num_io_procs set to 0')
      num_io_procs = 0
    END IF

    IF (lhave_radar) THEN
      IF (num_io_procs_radar /= 0) THEN
        CALL print_info_stderr(method_name, &
             & 'num_io_procs_radar has no effect if the model is compiled with the NOMPI compiler directive')
        CALL print_info_stderr(method_name, &
             & '--> num_io_procs_radar set to 0')
        num_io_procs_radar = 0
      ENDIF
    ENDIF

    IF (num_restart_procs /= 0) THEN
      CALL print_info_stderr(method_name, &
       & 'num_restart_procs has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
       & '--> num_restart_procs set to 0')
      num_restart_procs = 0
    END IF
    IF (sizeof_prefetch_processes /= 0) THEN
      CALL print_info_stderr(method_name, &
      & 'sizeof_prefetch_processes has no effect if the model is compiled with the NOMPI compiler directive')
      CALL print_info_stderr(method_name, &
      & '--> sizeof_prefetch_processes set to 0')
      sizeof_prefetch_processes = 0
    END IF


    ! set the sequential values
    p_work_pe0 = 0
    num_io_procs = 0
    sizeof_prefetch_processes = 0
    num_work_procs = 1
    p_comm_io = mpi_comm_self

    num_radario_procs = 0
    p_radario_pe0 = 0


    ! fill some derived variables
    process_mpi_all_workroot_id     = p_work_pe0
    process_mpi_all_ioroot_id       = p_io_pe0
    process_mpi_all_test_id         = p_work_pe0-1
    process_mpi_restart_size        = num_restart_procs
    process_mpi_io_size             = num_io_procs
    process_mpi_radario_size        = num_radario_procs
    process_mpi_all_radarioroot_id  = p_radario_pe0
    process_mpi_all_prefroot_id     = p_pref_pe0
    process_mpi_pref_size           = sizeof_prefetch_processes
    p_comm_work_pref_compute_pe0    = p_work_pe0

    ! In case of test run, only the test process is stdio
    process_is_stdio = (my_process_mpi_all_id == process_mpi_stdio_id)
    process_is_mpi_parallel = p_n_work > 1

    p_work_root_processes(1)%comp_id = comp_id
    p_work_root_processes(1)%global_mpi_id = my_global_mpi_id

! if there is a hamocc proccess, create hamocc-ocean communicator
    

    ! still to be filled
!     process_mpi_local_comm  = process_mpi_all_comm
!     process_mpi_local_size  = process_mpi_all_size
!     my_process_mpi_local_id = my_process_mpi_all_id

   


  END SUBROUTINE set_mpi_work_communicators
  !-------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  INTEGER FUNCTION get_mpi_work_intercomm(other_comp_id)

    INTEGER, INTENT(IN) :: other_comp_id
    CHARACTER(len=*), PARAMETER :: method_name = 'get_mpi_work_intercomm'

    INTEGER :: peer_comm, num_component, other_comp_root_global_mpi_id, i, &
               intercomm

    get_mpi_work_intercomm = -1

  END FUNCTION get_mpi_work_intercomm
  !-------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  SUBROUTINE set_default_mpi_work_variables()

    ! fill some derived variables
    process_is_stdio = (my_process_mpi_all_id == process_mpi_stdio_id)
    process_is_mpi_parallel = (process_mpi_all_size > 1)
    my_mpi_function = work_mpi_process

    process_mpi_all_test_id     = -1
    process_mpi_all_workroot_id = 0
    process_mpi_io_size         = 0
    process_mpi_restart_size    = 0
    process_mpi_pref_size       = 0
    process_mpi_radario_size    = 0
    p_comm_work_pref_compute_pe0 = 0

!     process_mpi_local_comm  = process_mpi_all_comm
!     process_mpi_local_size  = process_mpi_all_size
!     my_process_mpi_local_id = my_process_mpi_all_id

    ! set some of the old variables
    ! should be removed once the old variables are cleaned
    p_pe           = my_process_mpi_all_id
    p_io           = 0
    num_test_procs = 0
    num_work_procs = process_mpi_all_size
    p_work_pe0     = 0
    p_io_pe0       = process_mpi_all_size    ! Number of I/O PE 0 within all PEs (process_mpi_all_size if no I/O PEs)
    p_radario_pe0  = process_mpi_all_size
    ! Number of prefetching PE 0 within all PEs (process_mpi_all_size if no prefetching PEs)
    p_pref_pe0     = process_mpi_all_size
    p_n_work       = process_mpi_all_size
    p_pe_work      = my_process_mpi_all_id

    p_comm_work             = process_mpi_all_comm
    p_comm_work_io          = MPI_COMM_NULL
    p_comm_work_test        = MPI_COMM_NULL
    p_comm_work_2_io        = MPI_COMM_NULL
    p_comm_work_restart     = MPI_COMM_NULL
    p_comm_work_2_restart   = MPI_COMM_NULL
    p_comm_io               = MPI_COMM_NULL
    P_comm_work_pref        = MPI_COMM_NULL
    P_comm_work_2_pref      = MPI_COMM_NULL

    ! print some info
    IF ( .NOT. process_is_mpi_parallel) THEN
      WRITE (nerr,'(/,a,/)')  ' Single processor run.'
    ELSEIF (process_is_stdio) THEN
      WRITE (nerr,'(/,a,a,a,i0,a,/)') ' ', TRIM(process_mpi_name), &
        ' runs on ', process_mpi_all_size, ' mpi processes.'
    END IF

  END SUBROUTINE set_default_mpi_work_variables
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! Splits the global communicator into this component's communicator
  !! Should be called before the component configuration
  SUBROUTINE split_global_mpi_communicator(component_no, num_components)
    INTEGER, INTENT(in) :: component_no, num_components

    INTEGER :: new_communicator
    LOGICAL             :: l_mpi_is_initialised
    CHARACTER(len=*), PARAMETER :: &
         routine = modname//'::split_process_mpi_communicator'
    ALLOCATE(p_work_root_processes(1))
    RETURN

  END SUBROUTINE split_global_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  !>
  !! Set this component's communicator
  !! The new communicator id is duplicated from the input communicator
  !! Should be called before the component configuration
  SUBROUTINE set_process_mpi_communicator(new_communicator)
    INTEGER, INTENT(in) :: new_communicator

    LOGICAL             :: l_mpi_is_initialised
    CHARACTER(len=*), PARAMETER :: &
         routine = modname//'::set_process_mpi_communicator'

    process_mpi_all_comm    = new_communicator
    process_mpi_all_size    = 1
    my_process_mpi_all_id   = 0


    CALL set_default_mpi_work_variables()

  END SUBROUTINE set_process_mpi_communicator
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE start_mpi(global_name)

!#ifdef _OPENMP
!    USE mo_util_string, ONLY: toupper
!#endif

    CHARACTER(len=*), INTENT(in), OPTIONAL :: global_name

    ! variables are required for determing I/O size in bytes of the defined
    ! KIND types for assigning the right MPI data types with the used kinds
    INTEGER     :: iig = 0
    INTEGER(i4) :: ii4 = 0_i4
    INTEGER(i8) :: ii8 = 0_i8
    REAL        :: rrg = 0.0
    REAL(sp)    :: rsp = 0.0_sp
    REAL(dp)    :: rdp = 0.0_dp
    INTEGER :: version, subversion

    CHARACTER(len=132) :: yname

    CHARACTER(len=256) :: group_names(3) = ["icon ", "yac  ", "comin"]
    INTEGER            :: group_comms(3)

    ! variables used for determing the OpenMP threads
    ! suitable as well for coupled models

    CHARACTER(len=*), PARAMETER :: routine = modname//'::start_mpi'


! #ifdef _OPENMP
!     global_no_of_threads = 1
! #endif


    ! start MPI

    WRITE (nerr,'(a,a)')  routine, ' No MPI: Single processor run.'
    ! set defaults for sequential run
    global_mpi_size  = 1        ! total number of processes in global world
    my_global_mpi_id = 0        ! process id in global world
    is_global_mpi_parallel = .false.
    process_mpi_name = 'uknown'
    process_mpi_all_comm = MPI_COMM_NULL
    p_mpi_comm_null = MPI_COMM_NULL



    is_global_mpi_parallel = global_mpi_size > 1



    ! The number of threads, if varying, will be defined via
    ! namelists
! #ifdef _OPENMP
!     ! Expect that PE 0 did got the information of OMP_NUM_THREADS.
!     ! That might be wrong in the coupled case when the model is
!     ! started via MPI dynamic process creation. So we have to check
!     ! the environment variable too.
!     IF (my_global_mpi_id == 0) THEN
!
!       IF (is_global_mpi_parallel) THEN
!         WRITE (nerr,'(/,a)') ' Running globally hybrid OpenMP-MPI mode.'
!       ELSE
!         WRITE (nerr,'(/,a)') ' Running globally OpenMP mode.'
!       ENDIF
!
!       env_name = toupper(TRIM(yname)) // '_THREADS'
! #ifdef __SX__
!       CALL getenv(TRIM(env_name), thread_num)
!
!       IF (thread_num /= ' ') THEN
! #else
!       CALL get_environment_variable(name=TRIM(env_name), value=thread_num, &
!           status=istat)
!       IF (istat == 0) THEN
! #endif
!          READ(thread_num,*) global_no_of_threads
!        ELSE
!          WRITE (nerr,'(1x,a,/)') ' Global number of OpenMP threads not given!'
!          WRITE (nerr,'(1x,a,a,a,/,1x,a,a,a)') &
!              ' Environment variable ', TRIM(env_name), ' either not set,', &
!              ' or not available to ', TRIM(yname), ' root PE.'
!          global_no_of_threads = omp_get_max_threads()
!        ENDIF
!     ENDIF ! (my_global_mpi_id == 0)
!
! #ifndef 1
!     ! Make number of threads from environment available to all model PEs
!     CALL MPI_BCAST (global_no_of_threads, 1, MPI_INTEGER, 0, global_mpi_communicator, p_error)
! #endif
!     ! Inform on OpenMP thread usage
! !     CALL OMP_SET_NUM_THREADS(threads)
! !     threads = OMP_GET_MAX_THREADS()
!
!     IF (my_global_mpi_id == 0) THEN
!        WRITE (nerr,*)
!        WRITE (nerr,'(1x,a,i3)') ' global_no_of_threads is ', global_no_of_threads
!     ENDIF
! #endif


    ! by default, the global communicator is the process communicator
    CALL set_process_mpi_name(global_mpi_name)
    CALL set_process_mpi_communicator(global_mpi_communicator)

    CALL init_logger(proc_id=get_my_global_mpi_id(), &
                     l_write_output=my_process_is_stdio(), &
                     nerr_unit=nerr, &
                     l_extra_output=(proc_split .AND. comm_lev > 0 .AND. get_glob_proc0() == p_pe_work),&
                     extra_info_prefix='PROC SPLIT', &
                     callback_abort=abort_mpi)

  END SUBROUTINE start_mpi
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE stop_mpi

    INTEGER :: iexit = 0    
    ! finish MPI and clean up all PEs

    CALL finalize_coupler()

    CALL util_exit(iexit)
  END SUBROUTINE stop_mpi
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  SUBROUTINE abort_mpi

    ! this routine should be used instead of abort, util_abort() or
    ! STOP or any other exit call in all routines for proper clean up
    ! of all PEs


     CALL util_exit(1)

  END SUBROUTINE abort_mpi
  !------------------------------------------------------------------------------

  ! communicator set up
  SUBROUTINE p_set_communicator (nproca, nprocb, mapmesh, debug_parallel)

    INTEGER, INTENT(in) :: nproca, nprocb
    INTEGER, INTENT(in) :: mapmesh(0:,0:)
    INTEGER, INTENT(in) :: debug_parallel

  END SUBROUTINE p_set_communicator

!=========================================================================

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Comm_size()
  !---------------------------------------------------------------------------------------------------------------------------------
  FUNCTION p_comm_size(communicator)
    INTEGER :: p_comm_size
    INTEGER, INTENT(IN) :: communicator

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::p_comm_size'
    INTEGER :: ierr

    p_comm_size = 1
  END FUNCTION p_comm_size

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Comm_rank()
  !---------------------------------------------------------------------------------------------------------------------------------
  FUNCTION p_comm_rank(communicator)
    INTEGER :: p_comm_rank
    INTEGER, INTENT(IN) :: communicator

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::p_comm_rank'
    INTEGER :: ierr

    p_comm_rank = 0
  END FUNCTION p_comm_rank


!=========================================================================

  ! send implementation
  SUBROUTINE p_send_real (t_buffer, p_destination, p_tag, p_count, comm, use_g2g)

    REAL (dp), INTENT(in) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_send_real


  ! send implementation

  SUBROUTINE p_send_sreal (t_buffer, p_destination, p_tag, p_count, comm, use_g2g)

    REAL (sp), INTENT(in) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_send_sreal


  SUBROUTINE p_send_real_1d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_real_1d

  SUBROUTINE p_send_sreal_1d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (sp), INTENT(in) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_sreal_1d

  SUBROUTINE p_send_real_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_real_2d

  SUBROUTINE p_send_real_3d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_real_3d

  SUBROUTINE p_send_real_4d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_real_4d

  SUBROUTINE p_send_real_5d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(in) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_real_5d

  SUBROUTINE p_send_int (t_buffer, p_destination, p_tag, p_count, comm, use_g2g)

    INTEGER, INTENT(in) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_send_int

  SUBROUTINE p_send_int_1d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_int_1d

  SUBROUTINE p_send_int_2d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_int_2d

  SUBROUTINE p_send_int_3d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_int_3d

  SUBROUTINE p_send_int_4d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_int_4d


  SUBROUTINE p_send_bool (t_buffer, p_destination, p_tag, p_count, comm, use_g2g)

    LOGICAL, INTENT(in) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_send_bool

  SUBROUTINE p_send_bool_1d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_bool_1d


  SUBROUTINE p_send_char_1d (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(LEN=*), INTENT(in) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_char_1d


  SUBROUTINE p_send_bool_2d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_bool_2d

  SUBROUTINE p_send_bool_3d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_bool_3d

  SUBROUTINE p_send_bool_4d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(in) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_bool_4d

  SUBROUTINE p_send_char (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(in) :: t_buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_send_char

! non-blocking sends

  SUBROUTINE p_inc_request
    INTEGER, ALLOCATABLE :: tmp(:)


  END SUBROUTINE p_inc_request


  SUBROUTINE p_isend_real (t_buffer, p_destination, p_tag, p_count, comm, request, use_g2g)

    REAL (dp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_isend_real


  SUBROUTINE p_isend_sreal (t_buffer, p_destination, p_tag, p_count, comm, request, use_g2g)

    REAL (sp), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_isend_sreal


  SUBROUTINE p_isend_real_1d (t_buffer, p_destination, p_tag, p_count, comm, request)

    REAL (dp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request

  END SUBROUTINE p_isend_real_1d


  SUBROUTINE p_isend_sreal_1d (t_buffer, p_destination, p_tag, p_count, comm, request)

    REAL (sp), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request

  END SUBROUTINE p_isend_sreal_1d


  SUBROUTINE p_isend_real_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_real_2d


  SUBROUTINE p_isend_sreal_2d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (sp), INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_sreal_2d


  SUBROUTINE p_isend_real_3d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_real_3d

  SUBROUTINE p_isend_real_4d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_real_4d

  SUBROUTINE p_isend_real_5d (t_buffer, p_destination, p_tag, p_count, comm)

    REAL (dp), INTENT(inout) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_real_5d

  SUBROUTINE p_isend_int(t_buffer, p_destination, p_tag, p_count, comm, request, use_g2g)

    INTEGER, INTENT(in) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g


  END SUBROUTINE p_isend_int

  SUBROUTINE p_isend_int_1d(t_buffer, p_destination, p_tag, p_count, comm, request)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER,   INTENT(INOUT), OPTIONAL :: request


  END SUBROUTINE p_isend_int_1d

  SUBROUTINE p_isend_int_2d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_int_2d

  SUBROUTINE p_isend_int_3d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_isend_int_3d

  SUBROUTINE p_isend_int_4d (t_buffer, p_destination, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_int_4d


  SUBROUTINE p_isend_bool (t_buffer, p_destination, p_tag, p_count, comm, use_g2g)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g


  END SUBROUTINE p_isend_bool

  SUBROUTINE p_isend_bool_1d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_bool_1d

  SUBROUTINE p_isend_bool_2d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_bool_2d

  SUBROUTINE p_isend_bool_3d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_bool_3d

  SUBROUTINE p_isend_bool_4d (t_buffer, p_destination, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_bool_4d

  SUBROUTINE p_isend_char (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(inout) :: t_buffer
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_isend_char

  ! recv implementation

  SUBROUTINE p_recv_real (t_buffer, p_source, p_tag, p_count, comm, use_g2g)

    REAL (dp), INTENT(out) :: t_buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_recv_real

  ! recv implementation

  SUBROUTINE p_recv_sreal (t_buffer, p_source, p_tag, p_count, comm, use_g2g)

    REAL (sp), INTENT(out) :: t_buffer
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_recv_sreal

  SUBROUTINE p_isend_char_1d(t_buffer, p_destination, p_tag, p_count, comm)
    CHARACTER(len=*),  INTENT(inout) :: t_buffer(:)
    INTEGER,           INTENT(in) :: p_destination, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_isend_char_1d

  SUBROUTINE p_recv_real_1d (t_buffer, p_source, p_tag, p_count, comm, displs)

    REAL (dp), INTENT(out) :: t_buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm, displs

  END SUBROUTINE p_recv_real_1d

  SUBROUTINE p_recv_sreal_1d (t_buffer, p_source, p_tag, p_count, comm, displs)

    REAL (sp), INTENT(out) :: t_buffer(:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm, displs

  END SUBROUTINE p_recv_sreal_1d

  SUBROUTINE p_recv_real_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_real_2d

  SUBROUTINE p_recv_real_3d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_real_3d

  SUBROUTINE p_recv_real_4d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_real_4d

  SUBROUTINE p_recv_real_5d (t_buffer, p_source, p_tag, p_count, comm)

    REAL (dp), INTENT(out) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_real_5d

  SUBROUTINE p_recv_int (t_buffer, p_source, p_tag, p_count, comm, use_g2g)

    INTEGER, INTENT(out) :: t_buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_recv_int

  SUBROUTINE p_recv_int_1d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_int_1d

  SUBROUTINE p_recv_int_2d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_int_2d

  SUBROUTINE p_recv_int_3d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_int_3d

  SUBROUTINE p_recv_int_4d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_int_4d


  SUBROUTINE p_recv_bool (t_buffer, p_source, p_tag, p_count, comm, use_g2g)

    LOGICAL, INTENT(out) :: t_buffer
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_recv_bool

  SUBROUTINE p_recv_bool_1d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_bool_1d


  SUBROUTINE p_recv_char_1d (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(LEN=*), INTENT(out) :: t_buffer(:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_char_1d


  SUBROUTINE p_recv_bool_2d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_bool_2d

  SUBROUTINE p_recv_bool_3d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_bool_3d

  SUBROUTINE p_recv_bool_4d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(out) :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_bool_4d

  SUBROUTINE p_recv_char (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*),  INTENT(out) :: t_buffer
    INTEGER,           INTENT(in)  :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_recv_char

  ! non-blocking receives

  !================================================================================================
  ! CHARACTER SECTION -----------------------------------------------------------------------------
  !
  SUBROUTINE p_irecv_char (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER(len=*), INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_irecv_char

  SUBROUTINE p_irecv_char_1d(t_buffer, p_source, p_tag, p_count, comm, request)
    CHARACTER(len=*), INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER, OPTIONAL, INTENT(out) :: request
  END SUBROUTINE p_irecv_char_1d

  !================================================================================================
  ! REAL SECTION ----------------------------------------------------------------------------------
  !

  SUBROUTINE p_irecv_real (t_buffer, p_source, p_tag, p_count, comm, use_g2g)

    REAL(dp),  INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_irecv_real


  SUBROUTINE p_irecv_sreal (t_buffer, p_source, p_tag, p_count, comm, use_g2g)

    REAL(sp),  INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_irecv_sreal

  SUBROUTINE p_irecv_real_1d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_irecv_real_1d


  SUBROUTINE p_irecv_sreal_1d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(sp),  INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_irecv_sreal_1d


  SUBROUTINE p_irecv_real_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_real_2d


  SUBROUTINE p_irecv_sreal_2d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(sp),  INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_sreal_2d


  SUBROUTINE p_irecv_real_3d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_real_3d

  SUBROUTINE p_irecv_real_4d (t_buffer, p_source, p_tag, p_count, comm)

    REAL(dp),  INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_real_4d
  !================================================================================================
  ! INTEGER SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE p_irecv_int (t_buffer, p_source, p_tag, p_count, comm, request, use_g2g)

    INTEGER, INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER, OPTIONAL, INTENT(INOUT) :: request
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_irecv_int

  SUBROUTINE p_irecv_int_1d(t_buffer, p_source, p_tag, p_count, comm, &
    &                       request)

    INTEGER, INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    INTEGER, OPTIONAL, INTENT(out) :: request

  END SUBROUTINE p_irecv_int_1d

  SUBROUTINE p_irecv_int_2d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_int_2d

  SUBROUTINE p_irecv_int_3d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_int_3d

  SUBROUTINE p_irecv_int_4d (t_buffer, p_source, p_tag, p_count, comm)

    INTEGER, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_int_4d

  !================================================================================================
  ! LOGICAL SECTION -------------------------------------------------------------------------------
  !
  SUBROUTINE p_irecv_bool (t_buffer, p_source, p_tag, p_count, comm, use_g2g)

    LOGICAL, INTENT(inout) :: t_buffer
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

  END SUBROUTINE p_irecv_bool

  SUBROUTINE p_irecv_bool_1d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm

  END SUBROUTINE p_irecv_bool_1d

  SUBROUTINE p_irecv_bool_2d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_bool_2d

  SUBROUTINE p_irecv_bool_3d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_bool_3d

  SUBROUTINE p_irecv_bool_4d (t_buffer, p_source, p_tag, p_count, comm)

    LOGICAL, INTENT(inout) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in) :: p_source, p_tag
    INTEGER, OPTIONAL, INTENT(in) :: p_count, comm


  END SUBROUTINE p_irecv_bool_4d

  SUBROUTINE p_pack_int (t_var, t_buffer, p_pos, comm)

    INTEGER,   INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_pack_int

  SUBROUTINE p_pack_bool (t_var, t_buffer, p_pos, comm)

    LOGICAL,   INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_pack_bool

  SUBROUTINE p_pack_real (t_var, t_buffer, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_pack_real

  SUBROUTINE p_pack_int_1d (t_var, p_count, t_buffer, p_pos, comm)

    INTEGER,   INTENT(IN)    :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_pack_int_1d

  SUBROUTINE p_pack_real_1d (t_var, p_count, t_buffer, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_pack_real_1d


  SUBROUTINE p_pack_string (t_var, t_buffer, p_pos, comm)

    CHARACTER(LEN=*), INTENT(IN) :: t_var
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_pack_string

  SUBROUTINE p_pack_real_2d (t_var, p_count, t_buffer, p_pos, comm)

    REAL(wp),  INTENT(IN)    :: t_var(:,:)
    INTEGER,   INTENT(IN)    :: p_count
    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_pack_real_2d

  FUNCTION p_pack_size_int(p_count, comm) RESULT(pack_size)
    INTEGER, INTENT(in) :: p_count, comm
    INTEGER :: pack_size
    ! packing is only supported when mpi is available
    pack_size = -1
  END FUNCTION p_pack_size_int

  FUNCTION p_pack_size_bool(p_count, comm) result(pack_size)
    INTEGER, INTENT(in) :: p_count, comm
    INTEGER :: pack_size
    ! packing is only supported when mpi is available
    pack_size = -1
  END FUNCTION p_pack_size_bool

  FUNCTION p_pack_size_real_dp(p_count, comm) RESULT(pack_size)
    INTEGER, INTENT(in) :: p_count, comm
    INTEGER :: pack_size
    ! packing is only supported when mpi is available
    pack_size = -1
  END FUNCTION p_pack_size_real_dp

  FUNCTION p_pack_size_string(maxlen, comm) RESULT(pack_size)
    INTEGER, INTENT(in) :: maxlen, comm
    INTEGER :: pack_size, pack_size_int
    ! packing is only supported when mpi is available
    pack_size = -1
  END FUNCTION p_pack_size_string

  SUBROUTINE p_unpack_int (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER,   INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_unpack_int

  SUBROUTINE p_unpack_bool (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    LOGICAL,   INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_unpack_bool

  SUBROUTINE p_unpack_real (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(OUT)   :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_unpack_real

  SUBROUTINE p_unpack_int_1d (t_buffer, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    INTEGER,   INTENT(INOUT) :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_unpack_int_1d

  SUBROUTINE p_unpack_real_1d (t_buffer, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(INOUT) :: t_var(:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_unpack_real_1d


  SUBROUTINE p_unpack_string (t_buffer, p_pos, t_var, comm)

    CHARACTER, INTENT(IN) :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    CHARACTER(LEN=*), INTENT(INOUT) :: t_var
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_unpack_string


  SUBROUTINE p_unpack_real_2d (t_buffer, p_pos, t_var, p_count, comm)

    CHARACTER, INTENT(IN)    :: t_buffer(:)
    INTEGER,   INTENT(INOUT) :: p_pos
    REAL(wp),  INTENT(INOUT) :: t_var(:,:)
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
  END SUBROUTINE p_unpack_real_2d

  SUBROUTINE p_recv_packed (t_buffer, p_source, p_tag, p_count, comm)

    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_source, p_tag
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm

  END SUBROUTINE p_recv_packed

  SUBROUTINE p_irecv_packed (t_buffer, p_source, p_tag, p_count, comm, &
       request)

    CHARACTER, INTENT(INOUT) :: t_buffer(:)
    INTEGER,   INTENT(IN)    :: p_source, p_tag
    INTEGER,   INTENT(IN)    :: p_count
    INTEGER, OPTIONAL, INTENT(IN) :: comm
    INTEGER, OPTIONAL, INTENT(out) :: request

  END SUBROUTINE p_irecv_packed

  SUBROUTINE p_send_packed (t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER, INTENT(in) :: t_buffer(:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER,   INTENT(in) :: p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

  END SUBROUTINE p_send_packed

  SUBROUTINE p_send_packed_2d(t_buffer, p_destination, p_tag, p_count, comm)

    CHARACTER, INTENT(in) :: t_buffer(:,:)
    INTEGER,   INTENT(in) :: p_destination, p_tag
    INTEGER,   INTENT(in) :: p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm

  END SUBROUTINE p_send_packed_2d

  SUBROUTINE p_bcast_packed (t_buffer, p_source, p_count, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    CHARACTER :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER,   INTENT(in)    :: p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm
  END SUBROUTINE p_bcast_packed

  !
  !================================================================================================

  ! sendrecv implementation

  SUBROUTINE p_sendrecv_real_1d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm


  END SUBROUTINE p_sendrecv_real_1d

  SUBROUTINE p_sendrecv_real_2d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm


  END SUBROUTINE p_sendrecv_real_2d

  SUBROUTINE p_sendrecv_real_3d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm


  END SUBROUTINE p_sendrecv_real_3d

  SUBROUTINE p_sendrecv_real_4d (sendbuf, p_dest, recvbuf, p_source, &
                                  p_tag, comm)

    REAL(dp), INTENT(in)           :: sendbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_dest
    REAL(dp), INTENT(out)          :: recvbuf (:,:,:,:)
    INTEGER,  INTENT(in)           :: p_source
    INTEGER,  INTENT(in)           :: p_tag
    INTEGER,  INTENT(in) ,OPTIONAL :: comm


  END SUBROUTINE p_sendrecv_real_4d

  SUBROUTINE p_sendrecv_char_array (sendbuf, p_dest, recvbuf, p_source, p_tag, comm)
    CHARACTER(KIND = C_CHAR), INTENT(in) :: sendbuf (:)
    CHARACTER(KIND = C_CHAR), INTENT(out) :: recvbuf (:)
    INTEGER, VALUE :: p_dest, p_source, p_tag
    INTEGER, INTENT(in), OPTIONAL :: comm


  END SUBROUTINE p_sendrecv_char_array

  ! bcast implementation

  SUBROUTINE p_bcast_real (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(dp) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

  END SUBROUTINE p_bcast_real

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Bcast
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_bcast_real_single (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(sp) :: t_buffer
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

  END SUBROUTINE p_bcast_real_single

  SUBROUTINE p_bcast_real_1d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(dp) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

  END SUBROUTINE p_bcast_real_1d

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Bcast
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_bcast_real_1d_single (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(sp) :: t_buffer(:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm

  END SUBROUTINE p_bcast_real_1d_single

  SUBROUTINE p_bcast_real_2d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(dp) :: t_buffer(:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_real_2d


  SUBROUTINE p_bcast_real_2d_single (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(sp) :: t_buffer(:,:) ! SINGLE PRECISION
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_real_2d_single


  SUBROUTINE p_bcast_real_3d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL (dp) :: t_buffer(:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_real_3d

  SUBROUTINE p_bcast_real_4d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(dp) :: t_buffer(:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_real_4d

  SUBROUTINE p_bcast_real_5d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(dp) :: t_buffer(:,:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_real_5d

  SUBROUTINE p_bcast_real_7d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    REAL(dp) :: t_buffer(:,:,:,:,:,:,:)
    INTEGER,   INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_real_7d

  SUBROUTINE p_bcast_int_i4 (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER(i4) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_i4

  SUBROUTINE p_bcast_int_i8 (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER(i8) :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_i8

  SUBROUTINE p_bcast_int_1d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER :: t_buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_1d

  SUBROUTINE p_bcast_int_i8_1d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER(i8) :: t_buffer(:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_i8_1d

  SUBROUTINE p_bcast_int_2d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER :: t_buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_2d

  SUBROUTINE p_bcast_int_3d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER :: t_buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_3d

  SUBROUTINE p_bcast_int_4d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_4d

  SUBROUTINE p_bcast_int_7d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER                :: t_buffer(:,:,:,:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_int_7d

  SUBROUTINE p_bcast_bool (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    LOGICAL :: t_buffer
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_bool

  SUBROUTINE p_bcast_bool_1d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    LOGICAL :: t_buffer(:)
    INTEGER, INTENT(in) :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_bool_1d

  SUBROUTINE p_bcast_bool_2d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    LOGICAL :: t_buffer(:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_bool_2d

  SUBROUTINE p_bcast_bool_3d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    LOGICAL :: t_buffer(:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_bool_3d

  SUBROUTINE p_bcast_bool_4d (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    LOGICAL :: t_buffer(:,:,:,:)
    INTEGER, INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_bool_4d

  SUBROUTINE p_bcast_char (t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    CHARACTER(len=*) :: t_buffer
    INTEGER,           INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_bcast_char

  SUBROUTINE p_bcast_char_1d(t_buffer, p_source, comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    CHARACTER(*) :: t_buffer(:)
    INTEGER,       INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                      :: lexlength, flength


  END SUBROUTINE p_bcast_char_1d

  SUBROUTINE p_bcast_char_2d(t_buffer, p_source, comm)
    CHARACTER(*) :: t_buffer(:,:)
    INTEGER,       INTENT(in)    :: p_source
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER :: lexlength, flength

  END SUBROUTINE p_bcast_char_2d



  SUBROUTINE p_bcast_cchar (t_buffer, buflen, p_source, p_comm)
    ! intentionally without argument intent, because p_source determines if
    ! t_buffer is read or written
    INTEGER,           INTENT(IN)    :: buflen
    INTEGER(C_SIGNED_CHAR)                :: t_buffer(buflen)

    INTEGER,           INTENT(IN)    :: p_source
    INTEGER,           INTENT(IN)    :: p_comm


  END SUBROUTINE p_bcast_cchar

  ! A bcast variant that handles ALLOCATABLE strings. Cannot be overloaded with p_bcast() because that would make the CALL ambigous.
  ! After this CALL, the string will always be ALLOCATED, even IF its length IS zero
  !
  !XXX: Deactivated due to a gfortran bug that looses the contents of allocatable strings on return.
!  SUBROUTINE p_bcast_achar(string, source, comm)
!    CHARACTER(:), ALLOCATABLE, INTENT(INOUT) :: string
!    INTEGER, VALUE :: source, comm
!
!    INTEGER :: length, error
!    CHARACTER(*), PARAMETER :: routine = modname//":p_bcast_achar"
!
!#ifndef 1
!    ! inform the receivers about the length of the string
!    length = 0
!    IF(ALLOCATED(string)) length = LEN(string)
!    CALL p_bcast(length, source, comm)
!
!    ! ensure that the string IS ALLOCATED to the correct length
!    IF(ALLOCATED(string)) THEN
!        IF(length /= LEN(string) .OR. length == 0) DEALLOCATE(string)
!    END IF
!    IF(.NOT.ALLOCATED(string)) THEN
!        ALLOCATE(CHARACTER(LEN = length) :: string, STAT = error)
!        IF(error /= SUCCESS) CALL finish(routine, "memory allocation failed")
!    END IF
!
!    ! actually TRANSFER the string
!    CALL p_bcast(string, source, comm)
!#endif
!    END SUBROUTINE p_bcast_achar


!DR Test
  SUBROUTINE p_bcast_datetime (mtime_datetime, p_source, comm)

    TYPE(datetime), TARGET                  :: mtime_datetime
    INTEGER                 , INTENT(in)    :: p_source
    INTEGER       , OPTIONAL, INTENT(in)    :: comm
    !
    ! local
    character(len=max_datetime_str_len) :: mtime_datetime_str
    TYPE(datetime), POINTER             :: mtime_datetime_ptr
    TYPE(datetime), POINTER             :: datetime_loc
    INTEGER                             :: errno
 

  END SUBROUTINE p_bcast_datetime



  ! Collective CALL to determine whether this process IS a sender/receiver IN a broadcast operation.
  ! This routine IS robust IN the presence of inter-communicators (which IS its reason d'etre).
  SUBROUTINE p_get_bcast_role(root, comm, isSender, isReceiver)
    INTEGER, INTENT(IN) :: root, comm
    LOGICAL, INTENT(OUT) :: isSender, isReceiver
    isSender = .TRUE.
    isReceiver = .FALSE.
  END SUBROUTINE p_get_bcast_role

  ! probe implementation

  SUBROUTINE p_probe (p_tagcount, p_tagtable, p_source, &
       p_tag, p_count, comm)

    INTEGER,   INTENT(in)  :: p_tagcount, p_tagtable(:)
    INTEGER,   INTENT(out) :: p_source, p_tag, p_count
    INTEGER, OPTIONAL, INTENT(in) :: comm


  END SUBROUTINE p_probe
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE p_wait
  END SUBROUTINE p_wait

  SUBROUTINE p_wait_1(request)
    INTEGER, INTENT(INOUT) :: request
  END SUBROUTINE p_wait_1

  SUBROUTINE p_wait_n(requests)
    INTEGER, INTENT(INOUT) :: requests(:)
  END SUBROUTINE p_wait_n

  SUBROUTINE p_wait_n_stati(requests, stati)
    INTEGER, INTENT(INOUT) :: requests(:)
    INTEGER, INTENT(out) :: stati(:,:)
  END SUBROUTINE p_wait_n_stati

  SUBROUTINE p_wait_any(return_pe)

    INTEGER, INTENT(out) :: return_pe
    return_pe = 0
  END SUBROUTINE p_wait_any
  !------------------------------------------------------


  !------------------------------------------------------
  FUNCTION p_test() RESULT(ret)
    LOGICAL :: ret
    ret = .TRUE.
  END FUNCTION p_test


  !------------------------------------------------------
  SUBROUTINE p_barrier (comm)
  INTEGER ,INTENT(IN) ,OPTIONAL :: comm

  END SUBROUTINE p_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE global_mpi_barrier

  END SUBROUTINE global_mpi_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  SUBROUTINE work_mpi_barrier

  END SUBROUTINE work_mpi_barrier
  !------------------------------------------------------

  !------------------------------------------------------
  !> computes a global sum of real single precision numbers
  !
  !  perform an ALLREDUCE operation
  FUNCTION p_sum_sp_0d (zfield, comm) RESULT (p_sum)

    REAL(sp)                        :: p_sum
    REAL(sp),  INTENT(in)           :: zfield
    INTEGER,   INTENT(in)           :: comm
    p_sum = zfield
  END FUNCTION p_sum_sp_0d
  !------------------------------------------------------

  !------------------------------------------------------
  !> computes a global sum of real numbers
  !
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALLREDUCE operation
  FUNCTION p_sum_dp_0d (zfield, comm, root) RESULT (p_sum)

    REAL(dp)                        :: p_sum
    REAL(dp),  INTENT(in)           :: zfield
    INTEGER,   INTENT(in)           :: comm
    INTEGER,   INTENT(in), OPTIONAL :: root
    p_sum = zfield
  END FUNCTION p_sum_dp_0d
  !------------------------------------------------------

  !------------------------------------------------------
  FUNCTION p_sum_sp_1d (zfield, comm, root) RESULT (p_sum)

    REAL(sp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm, root
    REAL(sp)                      :: p_sum (SIZE(zfield))

    p_sum = zfield

  END FUNCTION p_sum_sp_1d

  !------------------------------------------------------
  FUNCTION p_sum_dp_1d (zfield, comm, root, use_g2g) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm, root
    REAL(dp)                      :: p_sum (SIZE(zfield))
    LOGICAL, OPTIONAL, INTENT(in) :: use_g2g
    LOGICAL :: loc_use_g2g

    p_sum = zfield

  END FUNCTION p_sum_dp_1d

  !------------------------------------------------------
  FUNCTION p_sum_dp_2d (zfield, comm, root) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm, root
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

    p_sum = zfield

  END FUNCTION p_sum_dp_2d

  !------------------------------------------------------
  FUNCTION p_sum_dp_3d (zfield, comm, root) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm, root
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2),SIZE(zfield,3))

    p_sum = zfield

  END FUNCTION p_sum_dp_3d

  FUNCTION p_sum_i8_1d (kfield, comm) RESULT (p_sum)

    INTEGER(i8),       INTENT(in) :: kfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER(i8)                   :: p_sum (SIZE(kfield))

    p_sum = kfield

  END FUNCTION p_sum_i8_1d

  FUNCTION p_sum_i_1d (kfield, comm, root) RESULT (p_sum)

    INTEGER,       INTENT(in) :: kfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm, root
    INTEGER                   :: p_sum (SIZE(kfield))

    p_sum = kfield

  END FUNCTION p_sum_i_1d


  FUNCTION p_sum_i_0d (kfield, comm) RESULT (p_sum)

    INTEGER,       INTENT(in) :: kfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                   :: p_sum

    p_sum = kfield

  END FUNCTION p_sum_i_0d

  FUNCTION p_global_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum
    REAL(dp)                      :: pe_sums(SIZE(zfield))

    p_sum = SUM(zfield)

  END FUNCTION p_global_sum_1d

  FUNCTION p_field_sum_1d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield))

    p_sum = zfield

  END FUNCTION p_field_sum_1d

  FUNCTION p_field_sum_2d (zfield, comm) RESULT (p_sum)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

    p_sum = zfield

  END FUNCTION p_field_sum_2d

  !> common code of min/max reductions
  SUBROUTINE p_minmax_common(in_field, out_field, n, op, loc_op, &
       proc_id, keyval, comm, root)
    INTEGER, INTENT(in) :: n, op, loc_op
    REAL(dp), INTENT(in) :: in_field(n)
    REAL(dp), INTENT(out) :: out_field(n)

    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(n)
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(n)
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

    out_field = in_field
  END SUBROUTINE p_minmax_common

  SUBROUTINE p_minmax_common_sp(in_field, out_field, n, op, loc_op, &
       proc_id, keyval, comm, root)
    INTEGER, INTENT(in) :: n, op, loc_op
    REAL(sp), INTENT(in) :: in_field(n)
    REAL(sp), INTENT(out) :: out_field(n)

    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(n)
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(n)
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

    out_field = in_field
  END SUBROUTINE p_minmax_common_sp


  !> computes a global maximum of real numbers
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_max_0d (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(dp)                         :: p_max
    REAL(dp),          INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id
    INTEGER, OPTIONAL, INTENT(inout) :: keyval
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

    REAL(dp) :: temp_in(1), temp_out(1)
    INTEGER :: temp_keyval(1), temp_proc_id(1)
    temp_in(1) = zfield
    IF (PRESENT(proc_id) .AND. PRESENT(keyval)) THEN
      temp_keyval(1) = keyval; temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           proc_id=temp_proc_id, keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1); proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(proc_id)) THEN
      temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           proc_id=temp_proc_id, comm=comm, root=root)
      proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(keyval)) THEN
      temp_keyval(1) = keyval
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
         keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1)
    ELSE ! .not. present(keyval) .and. .not. present(proc_id)
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           comm=comm, root=root)
    END IF

    p_max = temp_out(1)

  END FUNCTION p_max_0d

  FUNCTION p_max_0d_sp (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(sp)                         :: p_max
    REAL(sp),          INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id
    INTEGER, OPTIONAL, INTENT(inout) :: keyval
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

    REAL(dp) :: temp_in(1), temp_out(1)
    INTEGER :: temp_keyval(1), temp_proc_id(1)
    temp_in(1) = zfield
    IF (PRESENT(proc_id) .AND. PRESENT(keyval)) THEN
      temp_keyval(1) = keyval; temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           proc_id=temp_proc_id, keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1); proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(proc_id)) THEN
      temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           proc_id=temp_proc_id, comm=comm, root=root)
      proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(keyval)) THEN
      temp_keyval(1) = keyval
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
         keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1)
    ELSE ! .not. present(keyval) .and. .not. present(proc_id)
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_max, mpi_maxloc, &
           comm=comm, root=root)
    END IF

    p_max = temp_out(1)

  END FUNCTION p_max_0d_sp

  FUNCTION p_max_int_0d (zfield, comm) RESULT (p_max)

    INTEGER                          :: p_max
    INTEGER,           INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(in)    :: comm
    p_max = zfield

  END FUNCTION p_max_int_0d

  !> computes a global maximum of a real number array
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_max_1d (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(dp),          INTENT(in)    :: zfield(:)
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm
    REAL(dp)                         :: p_max (SIZE(zfield))

    CALL p_minmax_common(zfield, p_max, SIZE(zfield), mpi_max, mpi_maxloc, &
           proc_id=proc_id, keyval=keyval, comm=comm, root=root)

  END FUNCTION p_max_1d

  FUNCTION p_max_1d_sp (zfield, proc_id, keyval, comm, root) RESULT (p_max)

    REAL(sp),          INTENT(in)    :: zfield(:)
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm
    REAL(sp)                         :: p_max (SIZE(zfield))

    CALL p_minmax_common(zfield, p_max, SIZE(zfield), mpi_max, mpi_maxloc, &
           proc_id=proc_id, keyval=keyval, comm=comm, root=root)

  END FUNCTION p_max_1d_sp

  ! Computes maximum of a 1D field of integers.
  !
  ! @param[in] root Optional root PE, otherwise we perform an
  !                 ALL-TO-ALL operation
  FUNCTION p_max_int_1d (zfield, comm, root) RESULT (p_max)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER, OPTIONAL, INTENT(in) :: root
    INTEGER                       :: p_max (SIZE(zfield))

    p_max = zfield

  END FUNCTION p_max_int_1d

  FUNCTION p_max_2d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

    p_max = zfield

  END FUNCTION p_max_2d

  FUNCTION p_max_2d_sp (zfield, comm) RESULT (p_max)

    REAL(sp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(sp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2))

    p_max = zfield

  END FUNCTION p_max_2d_sp

  FUNCTION p_max_3d (zfield, comm) RESULT (p_max)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

    p_max = zfield

  END FUNCTION p_max_3d

  FUNCTION p_max_3d_sp (zfield, comm) RESULT (p_max)

    REAL(sp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(sp)                      :: p_max (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))

    p_max = zfield

  END FUNCTION p_max_3d_sp


  !> computes a global minimum of real numbers
  !
  ! @param[out]   proc_id  (Optional:) PE number of maximum value
  ! @param[inout] keyval   (Optional:) additional meta information
  ! @param[in]    root     (Optional:) root PE, otherwise we perform an
  !                                    ALL-TO-ALL operation
  !
  ! The parameter @p keyval can be used to communicate
  ! additional data on the maximum value, e.g., the level
  ! index where the maximum occurred.
  !
  FUNCTION p_min_0d (zfield, proc_id, keyval, comm, root) RESULT (p_min)

    REAL(dp)                         :: p_min
    REAL(dp),          INTENT(in)    :: zfield
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id
    INTEGER, OPTIONAL, INTENT(inout) :: keyval
    INTEGER, OPTIONAL, INTENT(in)    :: root
    INTEGER, OPTIONAL, INTENT(in)    :: comm

    REAL(dp) :: temp_in(1), temp_out(1)
    INTEGER :: temp_keyval(1), temp_proc_id(1)
    temp_in(1) = zfield
    IF (PRESENT(proc_id) .AND. PRESENT(keyval)) THEN
      temp_keyval(1) = keyval; temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
           proc_id=temp_proc_id, keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1); proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(proc_id)) THEN
      temp_proc_id(1) = proc_id
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
           proc_id=temp_proc_id, comm=comm, root=root)
      proc_id = temp_proc_id(1)
    ELSE IF (PRESENT(keyval)) THEN
      temp_keyval(1) = keyval
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
         keyval=temp_keyval, comm=comm, root=root)
      keyval = temp_keyval(1)
    ELSE ! .not. present(keyval) .and. .not. present(proc_id)
      CALL p_minmax_common(temp_in, temp_out, 1, mpi_min, mpi_minloc, &
           comm=comm, root=root)
    END IF

    p_min = temp_out(1)

  END FUNCTION p_min_0d

  FUNCTION p_min_int_0d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min
    p_min = zfield

  END FUNCTION p_min_int_0d

  FUNCTION p_min_1d (zfield, proc_id, keyval, comm, root) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(inout) :: proc_id(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(inout) :: keyval(SIZE(zfield))
    INTEGER, OPTIONAL, INTENT(in) :: root
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield))

    CALL p_minmax_common(zfield, p_min, SIZE(zfield), mpi_min, mpi_minloc, &
           proc_id=proc_id, keyval=keyval, comm=comm, root=root)

  END FUNCTION p_min_1d

  FUNCTION p_min_int_1d (zfield, comm) RESULT (p_min)

    INTEGER,           INTENT(in) :: zfield(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    INTEGER                       :: p_min (SIZE(zfield))

    p_min = zfield

  END FUNCTION p_min_int_1d

  FUNCTION p_min_2d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2))

    p_min = zfield

  END FUNCTION p_min_2d

  FUNCTION p_min_3d (zfield, comm) RESULT (p_min)

    REAL(dp),          INTENT(in) :: zfield(:,:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_min (SIZE(zfield,1),SIZE(zfield,2)&
                                           ,SIZE(zfield,3))
    p_min = zfield

  END FUNCTION p_min_3d

  FUNCTION p_lor_0d(zfield, comm) RESULT(res)
    LOGICAL :: res
    LOGICAL, INTENT(in) :: zfield
    INTEGER, OPTIONAL, INTENT(in) :: comm
    res = zfield
  END FUNCTION p_lor_0d

  ! Computes maximum of a 1D field of integers.
  !
  ! @param[in] root Optional root PE, otherwise we perform an
  !                 ALL-TO-ALL operation
  SUBROUTINE p_allreduce_max_int_1d (zfield, comm, root)
    INTEGER, INTENT(inout) :: zfield(:)
    INTEGER, INTENT(in) :: comm
    INTEGER, INTENT(in), OPTIONAL :: root

  END SUBROUTINE p_allreduce_max_int_1d

  FUNCTION p_reduce_i8_0d(send, op, root, comm) RESULT(recv)
    INTEGER(i8) :: recv
    INTEGER(i8), INTENT(IN) :: send
    INTEGER, INTENT(in) :: op, root, comm

    recv = send
  END FUNCTION p_reduce_i8_0d

  FUNCTION p_reduce_i4_0d(send, op, root, comm) &
       RESULT(recv)
    INTEGER(i4) :: recv
    INTEGER(i4), INTENT(IN) :: send
    INTEGER, INTENT(in) :: op, root, comm

    recv = send
  END FUNCTION p_reduce_i4_0d

  INTEGER(i4) FUNCTION p_allreduce_int4_0d(input, reductionOp, comm) RESULT(resultVar)
    INTEGER(i4), INTENT(IN) :: input
    INTEGER, INTENT(in) :: reductionOp, comm

    resultVar = input
  END FUNCTION p_allreduce_int4_0d

  FUNCTION p_allreduce_bool_0d(input, reductionOp, comm) &
       RESULT(res)
    LOGICAL :: res
    LOGICAL, INTENT(IN) :: input
    INTEGER, INTENT(in) :: reductionOp, comm

    res = input
  END FUNCTION p_allreduce_bool_0d

  INTEGER(i8) FUNCTION p_allreduce_int8_0d(input, reductionOp, comm) RESULT(resultVar)
    INTEGER(i8), INTENT(IN) :: input
    INTEGER, INTENT(in) :: reductionOp, comm

    resultVar = input
  END FUNCTION p_allreduce_int8_0d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatter
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_scatter_real_1d1d(sendbuf, recvbuf, p_src, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

     recvbuf = sendbuf
   END SUBROUTINE p_scatter_real_1d1d


   SUBROUTINE p_scatter_real_2d1d(sendbuf, recvbuf, p_src, comm)
    REAL(dp),          INTENT(inout) :: sendbuf(:,:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

    recvbuf = sendbuf(:,1)
  END SUBROUTINE p_scatter_real_2d1d

  SUBROUTINE p_scatter_sp_1d1d(sendbuf, recvbuf, p_src, comm)
    REAL(sp),          INTENT(inout) :: sendbuf(:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

     recvbuf = sendbuf
   END SUBROUTINE p_scatter_sp_1d1d

  SUBROUTINE p_scatter_sp_2d1d(sendbuf, recvbuf, p_src, comm)
    REAL(sp),          INTENT(inout) :: sendbuf(:,:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

    recvbuf = sendbuf(:,1)
  END SUBROUTINE p_scatter_sp_2d1d

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatter
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_scatter_int_1d1d(sendbuf, recvbuf, p_src, comm)
    INTEGER,           INTENT(inout) :: sendbuf(:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

     recvbuf = sendbuf
  END SUBROUTINE p_scatter_int_1d1d

  SUBROUTINE p_scatter_int_2d1d(sendbuf, recvbuf, p_src, comm)
    INTEGER,           INTENT(inout) :: sendbuf(:,:), recvbuf(:)
    INTEGER,           INTENT(in) :: p_src
    INTEGER, OPTIONAL, INTENT(in) :: comm

     recvbuf = sendbuf(:,1)
  END SUBROUTINE p_scatter_int_2d1d

  SUBROUTINE p_gather_real_0d1d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(in   ) :: sendbuf
    REAL(dp),          INTENT(inout) :: recvbuf(:)
    INTEGER,           INTENT(in   ) :: p_dest
    INTEGER, OPTIONAL, INTENT(in   ) :: comm

     recvbuf(:) = sendbuf
  END SUBROUTINE p_gather_real_0d1d

  SUBROUTINE p_gather_real_1d2d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(in   ) :: sendbuf(:)
    REAL(dp),          INTENT(inout) :: recvbuf(:,:)
    INTEGER,           INTENT(in   ) :: p_dest
    INTEGER, OPTIONAL, INTENT(in   ) :: comm

     recvbuf(:,1) = sendbuf(:)
  END SUBROUTINE p_gather_real_1d2d

  SUBROUTINE p_gather_real_2d3d(sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(in   ) :: sendbuf(:,:)
    REAL(dp),          INTENT(inout) :: recvbuf(:,:,:)
    INTEGER,           INTENT(in   ) :: p_dest
    INTEGER, OPTIONAL, INTENT(in   ) :: comm

    recvbuf(:,:,1) = sendbuf(:,:)
  END SUBROUTINE p_gather_real_2d3d

   SUBROUTINE p_gather_real_5d6d (sendbuf, recvbuf, p_dest, comm)
     REAL(dp),          INTENT(in   ) :: sendbuf(:,:,:,:,:)
     REAL(dp),          INTENT(inout) :: recvbuf(:,:,:,:,:,:)
     INTEGER,           INTENT(in   ) :: p_dest
     INTEGER, OPTIONAL, INTENT(in   ) :: comm

     recvbuf(:,:,:,:,:,LBOUND(recvbuf,6)) = sendbuf(:,:,:,:,:)
   END SUBROUTINE p_gather_real_5d6d


  SUBROUTINE p_gather_real_1d1d (sendbuf, recvbuf, p_dest, comm)
    REAL(dp),          INTENT(in   ) :: sendbuf(:)
    REAL(dp),          INTENT(inout) :: recvbuf(:)
    INTEGER,           INTENT(in   ) :: p_dest
    INTEGER, OPTIONAL, INTENT(in   ) :: comm

     recvbuf(:) = sendbuf(:)
   END SUBROUTINE p_gather_real_1d1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Gather()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_gather_int_0d1d (sendbuf, recvbuf, p_dest, comm)
     INTEGER,           INTENT(in   ) :: sendbuf
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in   ) :: p_dest
     INTEGER, OPTIONAL, INTENT(in   ) :: comm

     recvbuf = sendbuf
   END SUBROUTINE p_gather_int_0d1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Gather()
  !---------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE p_gather_int_1d1d (sendbuf, recvbuf, p_dest, comm)
     INTEGER,           INTENT(in   ) :: sendbuf(:)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in   ) :: p_dest
     INTEGER, OPTIONAL, INTENT(in   ) :: comm

     recvbuf = sendbuf
  END SUBROUTINE p_gather_int_1d1d

  SUBROUTINE p_gather_int_1d2d(sendbuf, recvbuf, p_dest, comm)
    INTEGER,           INTENT(inout) :: recvbuf(:,:)
    INTEGER,           INTENT(in) :: p_dest, sendbuf(:)
    INTEGER, OPTIONAL, INTENT(in) :: comm

     recvbuf(:,1) = sendbuf(:)
   END SUBROUTINE p_gather_int_1d2d

  SUBROUTINE p_gather_int_2d3d(sendbuf, recvbuf, p_dest, comm)
    INTEGER,           INTENT(in) :: sendbuf(:,:)
    INTEGER,        INTENT(inout) :: recvbuf(:,:,:)
    INTEGER,           INTENT(in) :: p_dest
    INTEGER, OPTIONAL, INTENT(in) :: comm

     recvbuf(:,:,1) = sendbuf(:,:)
   END SUBROUTINE p_gather_int_2d3d

  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Gather()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_gather_char_0d1d (sbuf, recvbuf, p_dest, comm)
     CHARACTER(len=*),  INTENT(in)    ::  sbuf
     CHARACTER(len=*),  INTENT(inout) ::  recvbuf(:)
     INTEGER,           INTENT(in) :: p_dest
     INTEGER, OPTIONAL, INTENT(in) :: comm
     CHARACTER(*), PARAMETER :: routine = modname//"::p_gather_char_0d1d"

     IF (LEN(sbuf) /= LEN(recvbuf(1))) THEN
       CALL finish (routine, 'Internal error: String lengths do not match!')
     END IF
     recvbuf = sbuf
   END SUBROUTINE p_gather_char_0d1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Gather()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_gather_bool_0d1d (sendbuf, recvbuf, p_dest, comm)
     LOGICAL,           INTENT(in   ) :: sendbuf
     LOGICAL,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in   ) :: p_dest
     INTEGER, OPTIONAL, INTENT(in   ) :: comm

     recvbuf = sendbuf
   END SUBROUTINE p_gather_bool_0d1d


   SUBROUTINE p_gatherv_int (sendbuf, sendcount, recvbuf, recvcounts, &
     &                       displs, p_dest, comm)
     INTEGER,           INTENT(in)    :: sendbuf(:), sendcount
     INTEGER,           INTENT(inout) :: recvbuf(:), recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER, OPTIONAL, INTENT(in)    :: comm

     recvbuf((displs(1)+1):(displs(1)+sendcount)) = sendbuf(1:sendcount)
   END SUBROUTINE p_gatherv_int


   SUBROUTINE p_gatherv_real2D2D (sendbuf, sendcount, recvbuf, recvcounts, &
     &                            displs, p_dest, comm)
     REAL(DP), INTENT(IN) :: sendbuf(:,:)
     INTEGER, INTENT(IN)  :: sendcount
     REAL(DP), INTENT(OUT) :: recvbuf(:,:)
     INTEGER, INTENT(IN)  :: recvcounts(:)
     INTEGER, INTENT(IN)  :: displs(:)
     INTEGER, INTENT(IN)  :: p_dest
     INTEGER, INTENT(IN)  :: comm

     recvbuf(:, (displs(1)+1):(displs(1)+sendcount)) = sendbuf(:, 1:sendcount)
   END SUBROUTINE p_gatherv_real2D2D


   SUBROUTINE p_gatherv_sreal2D2D (sendbuf, sendcount, recvbuf, recvcounts, &
     &                            displs, p_dest, comm)
     REAL(SP), INTENT(IN) :: sendbuf(:,:)
     INTEGER, INTENT(IN)  :: sendcount
     REAL(SP), INTENT(OUT) :: recvbuf(:,:)
     INTEGER, INTENT(IN)  :: recvcounts(:)
     INTEGER, INTENT(IN)  :: displs(:)
     INTEGER, INTENT(IN)  :: p_dest
     INTEGER, INTENT(IN)  :: comm

     recvbuf(:, (displs(1)+1):(displs(1)+sendcount)) = sendbuf(:, 1:sendcount)
   END SUBROUTINE p_gatherv_sreal2D2D


   SUBROUTINE p_gatherv_int2D2D (sendbuf, sendcount, recvbuf, recvcounts, &
     &                            displs, p_dest, comm)
     INTEGER, INTENT(IN)  :: sendbuf(:,:)
     INTEGER, INTENT(IN)  :: sendcount
     INTEGER, INTENT(OUT)  :: recvbuf(:,:)
     INTEGER, INTENT(IN)  :: recvcounts(:)
     INTEGER, INTENT(IN)  :: displs(:)
     INTEGER, INTENT(IN)  :: p_dest
     INTEGER, INTENT(IN)  :: comm

     recvbuf(:, (displs(1)+1):(displs(1)+sendcount)) = sendbuf(:, 1:sendcount)
   END SUBROUTINE p_gatherv_int2D2D


   SUBROUTINE p_gatherv_real2D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:,:)
     INTEGER,           INTENT(IN)    :: sendcount
     REAL(dp),          INTENT(INOUT) :: recvbuf(:)
     INTEGER,           intent(IN)    :: recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

     recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
   END SUBROUTINE p_gatherv_real2D1D


  SUBROUTINE p_gatherv_int2D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
    INTEGER,           INTENT(IN)    :: sendbuf(:,:)
    INTEGER,           INTENT(IN)    :: sendcount
    INTEGER,           INTENT(INOUT) :: recvbuf(:)
    INTEGER,           INTENT(IN)    :: recvcounts(:), displs(:)
    INTEGER,           INTENT(in)    :: p_dest
    INTEGER,           INTENT(in)    :: comm

    recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
  END SUBROUTINE p_gatherv_int2D1D


   SUBROUTINE p_gatherv_real3D1D (sendbuf, sendcount, recvbuf, recvcounts, displs, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:,:,:)
     INTEGER,           INTENT(IN)    :: sendcount
     REAL(dp),          INTENT(INOUT) :: recvbuf(:)
     INTEGER,           intent(IN)    :: recvcounts(:), displs(:)
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

     recvbuf(:) = RESHAPE(sendbuf, (/ SIZE(recvbuf) /) )
   END SUBROUTINE p_gatherv_real3D1D


   SUBROUTINE p_scatterv_real1D2D (sendbuf, sendcounts, displs, recvbuf, recvcount, p_dest, comm)
     REAL(dp),          INTENT(IN)    :: sendbuf(:)
     INTEGER,           INTENT(IN)    :: sendcounts(:), displs(:)
     REAL(dp),          INTENT(INOUT) :: recvbuf(:,:)
     INTEGER,           INTENT(IN)    :: recvcount
     INTEGER,           INTENT(in)    :: p_dest
     INTEGER,           INTENT(in)    :: comm

     recvbuf(:,:) = RESHAPE(sendbuf, (/ SIZE(recvbuf,1), SIZE(recvbuf,2) /))
   END SUBROUTINE p_scatterv_real1D2D


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatterv()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_scatterv_real1D1D (sendbuf, sendcounts, displs, recvbuf, recvcount, p_src, comm)
        implicit none
        REAL(wp), INTENT(IN) :: sendbuf(:)
        INTEGER, INTENT(IN)  :: sendcounts(:)
        REAL(wp), INTENT(INOUT) :: recvbuf(:)
        INTEGER, INTENT(IN)  :: recvcount
        INTEGER, INTENT(IN)  :: displs(:)
        INTEGER, INTENT(IN)  :: p_src
        INTEGER, INTENT(IN)  :: comm

        recvbuf(1:recvcount) = sendbuf((displs(1)+1):(displs(1)+recvcount))
   END SUBROUTINE p_scatterv_real1D1D

   SUBROUTINE p_scatterv_int_1d1d(sendbuf, sendcounts, displs, recvbuf, &
     recvcount, p_src, comm)
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: sendbuf(:)
     INTEGER, INTENT(IN) :: sendcounts(:)
     INTEGER, INTENT(INOUT) :: recvbuf(:)
     INTEGER, INTENT(IN)  :: recvcount
     INTEGER, INTENT(IN)  :: displs(:)
     INTEGER, INTENT(IN)  :: p_src
     INTEGER, INTENT(IN)  :: comm

     recvbuf(1:recvcount) = sendbuf((displs(1)+1):(displs(1)+recvcount))
   END SUBROUTINE p_scatterv_int_1d1d


  !---------------------------------------------------------------------------------------------------------------------------------
  !> wrapper for MPI_Scatterv()
  !---------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE p_scatterv_single1D1D (sendbuf, sendcounts, displs, recvbuf, recvcount, p_src, comm)
        implicit none
        REAL(sp), INTENT(IN) :: sendbuf(:)
        INTEGER, INTENT(IN)  :: sendcounts(:)
        REAL(sp), INTENT(INOUT) :: recvbuf(:)
        INTEGER, INTENT(IN)  :: recvcount
        INTEGER, INTENT(IN)  :: displs(:)
        INTEGER, INTENT(IN)  :: p_src
        INTEGER, INTENT(IN)  :: comm

        recvbuf(1:recvcount) = sendbuf((displs(1)+1):(displs(1)+recvcount))
   END SUBROUTINE p_scatterv_single1D1D

   SUBROUTINE p_allgather_int_0d1d(sendbuf, recvbuf, sendcount, recvcount, comm)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in) :: sendbuf
     INTEGER, OPTIONAL, INTENT(in) :: sendcount, recvcount, comm

     recvbuf = sendbuf
   END SUBROUTINE p_allgather_int_0d1d

   SUBROUTINE p_allgather_int_1d2d(sendbuf, recvbuf, sendcount, recvcount, comm)
     INTEGER,           INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: sendbuf(:)
     INTEGER, OPTIONAL, INTENT(in) :: sendcount, recvcount, comm

     recvbuf(:, 1) = sendbuf
   END SUBROUTINE p_allgather_int_1d2d

   SUBROUTINE p_allgatherv_real_1d(sendbuf, recvbuf, recvcounts, comm)
     REAL(dp),          INTENT(in)    :: sendbuf(:)
     REAL(dp),          INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in)    :: recvcounts(:)
     INTEGER, OPTIONAL, INTENT(in)    :: comm

     recvbuf = sendbuf
   END SUBROUTINE p_allgatherv_real_1d

   SUBROUTINE p_allgatherv_int_1d(sendbuf, recvbuf, recvcounts, displs, &
     &                            comm)
     INTEGER,           INTENT(in)    :: sendbuf(:)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in)    :: recvcounts(:), displs(:)
     INTEGER, OPTIONAL, INTENT(in)    :: comm

     recvbuf = sendbuf
   END SUBROUTINE p_allgatherv_int_1d

   SUBROUTINE p_allgatherv_int_1d_contiguous(sendbuf, recvbuf, recvcounts, &
     &                                       comm)
     INTEGER,           INTENT(in)    :: sendbuf(:)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in)    :: recvcounts(:)
     INTEGER, OPTIONAL, INTENT(in)    :: comm

     recvbuf = sendbuf
   END SUBROUTINE p_allgatherv_int_1d_contiguous


   ! Commits a user-defined MPI type
   !
   FUNCTION p_commit_type_struct(oldtypes, blockcounts) RESULT(newtype)
     INTEGER :: newtype
     INTEGER, INTENT(IN) :: oldtypes(2), blockcounts(2)
     newtype = 0
   END FUNCTION p_commit_type_struct


   SUBROUTINE p_alltoall_int (sendbuf, recvbuf, comm)
     INTEGER,           INTENT(inout) :: sendbuf(:), recvbuf(:)
     INTEGER,           INTENT(in) :: comm
     recvbuf(:) = sendbuf(:)
   END SUBROUTINE p_alltoall_int


   SUBROUTINE p_alltoallv_real_2d (sendbuf, sendcounts, sdispls, &
     &                             recvbuf, recvcounts, rdispls, comm)
     REAL(dp), TARGET,  INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     REAL(dp),          INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
   END SUBROUTINE p_alltoallv_real_2d


   SUBROUTINE p_alltoallv_sreal_2d (sendbuf, sendcounts, sdispls, &
     &                              recvbuf, recvcounts, rdispls, comm)
     REAL(sp), TARGET,  INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     REAL(sp),          INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
   END SUBROUTINE p_alltoallv_sreal_2d


   SUBROUTINE p_alltoallv_int_2d (sendbuf, sendcounts, sdispls, &
     &                            recvbuf, recvcounts, rdispls, comm)
     INTEGER, TARGET,   INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     INTEGER,           INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
   END SUBROUTINE p_alltoallv_int_2d


   SUBROUTINE p_alltoallv_int (sendbuf, sendcounts, sdispls, &
     &                         recvbuf, recvcounts, rdispls, comm)
     INTEGER,           INTENT(in) :: sendbuf(:), sendcounts(:), sdispls(:)
     INTEGER,           INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     ! displs are zero based -> have to add 1
     recvbuf(rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(sdispls(1)+1:sdispls(1)+sendcounts(1))
   END SUBROUTINE p_alltoallv_int

   SUBROUTINE p_alltoallv_int_i8_1d(sendbuf, sendcounts, sdispls, &
     &                         recvbuf, recvcounts, rdispls, comm)
     INTEGER(i8),       INTENT(in) :: sendbuf(:)
     INTEGER(i8),       INTENT(inout) :: recvbuf(:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:), &
       &                              recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     ! displs are zero based -> have to add 1
     recvbuf(rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(sdispls(1)+1:sdispls(1)+sendcounts(1))
   END SUBROUTINE p_alltoallv_int_i8_1d


   SUBROUTINE p_alltoallv_p2p_real_2d (sendbuf, sendcounts, sdispls, &
     &                             recvbuf, recvcounts, rdispls, comm)
     REAL(wp),          INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     REAL(wp),          INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
   END SUBROUTINE p_alltoallv_p2p_real_2d


   SUBROUTINE p_alltoallv_p2p_int_2d (sendbuf, sendcounts, sdispls, &
     &                                recvbuf, recvcounts, rdispls, comm)
     INTEGER,           INTENT(in) :: sendbuf(:,:)
     INTEGER,           INTENT(in) :: sendcounts(:), sdispls(:)
     INTEGER,           INTENT(inout) :: recvbuf(:,:)
     INTEGER,           INTENT(in) :: recvcounts(:), rdispls(:)
     INTEGER,           INTENT(in) :: comm
     ! displs are zero based -> have to add 1
     recvbuf(:,rdispls(1)+1:rdispls(1)+recvcounts(1)) = &
       sendbuf(:,sdispls(1)+1:sdispls(1)+sendcounts(1))
   END SUBROUTINE p_alltoallv_p2p_int_2d

  SUBROUTINE p_clear_request(request)
    INTEGER, INTENT(INOUT) :: request
    request = mpi_request_null
  END SUBROUTINE p_clear_request

  SUBROUTINE p_clear_requests(requests)
    INTEGER, INTENT(INOUT) :: requests(:)
    requests = mpi_request_null
  END SUBROUTINE p_clear_requests

  FUNCTION p_mpi_wtime()
    REAL(dp) :: p_mpi_wtime
    p_mpi_wtime =  0d0
  END FUNCTION p_mpi_wtime


  !--------------------------------------------------------------------


  !> @return Global MPI ranks within communicator "comm"
  !
  SUBROUTINE get_mpi_comm_world_ranks(comm, global_ranks, nranks)
    INTEGER, INTENT(IN)  :: comm               !< MPI communicator
    INTEGER, ALLOCATABLE, INTENT(OUT) :: global_ranks(:)    !< Output: list of global MPI ranks in communicator "comm"
    INTEGER, INTENT(OUT) :: nranks             !< Output: number of entries in rank list
    ! local variables
    CHARACTER(*), PARAMETER :: routine = modname//"::get_mpi_comm_world_ranks"
    INTEGER              :: p_error, grp_comm, grp_comm_world, i
    INTEGER, ALLOCATABLE :: comm_ranks(:)

    nranks = 1
    ALLOCATE(global_ranks(1))
    global_ranks(1) = 0
  END SUBROUTINE get_mpi_comm_world_ranks

  !--------------------------------------------------------------------

  LOGICAL FUNCTION p_comm_is_intercomm(intercomm)

    INTEGER, INTENT(IN) :: intercomm
    LOGICAL :: flag
    INTEGER :: p_error

    p_comm_is_intercomm = .FALSE.
  END FUNCTION p_comm_is_intercomm

  !--------------------------------------------------------------------

  INTEGER FUNCTION p_comm_remote_size(intercomm)

    INTEGER, INTENT(IN) :: intercomm
    INTEGER :: remote_size, p_error

    p_comm_remote_size = 0

  END FUNCTION p_comm_remote_size


  LOGICAL FUNCTION p_isEqual_int(val, comm) RESULT(resultVar)
    INTEGER, INTENT(IN) :: val
    INTEGER, OPTIONAL, INTENT(IN) :: comm

    resultVar = .TRUE.
  END FUNCTION p_isEqual_int

  LOGICAL FUNCTION p_isEqual_charArray(charArray, comm) RESULT(resultVar)
    CHARACTER(KIND = C_CHAR), INTENT(IN) :: charArray(:)
    INTEGER, OPTIONAL, INTENT(IN) :: comm

    resultVar = .TRUE.
  END FUNCTION p_isEqual_charArray

END MODULE mo_mpi
