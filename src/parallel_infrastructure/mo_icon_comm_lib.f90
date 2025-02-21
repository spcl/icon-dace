! A collection of MPI communication tools
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


!#include "dsl_definitions.inc"
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


!--------------------------------------------------
! timers definition
!needs:
!   USE mo_timer, ONLY: timer_start, timer_stop, timers_level, <timers_names>...
!



!----------------------------
MODULE mo_icon_comm_lib

  USE mo_kind,            ONLY: wp
  USE mo_exception,       ONLY: message_text, message, finish, warning
  USE mo_parallel_config, ONLY: icon_comm_debug, max_send_recv_buffer_size, &
    & icon_comm_method, icon_comm_openmp, max_no_of_comm_variables,         &
    & max_no_of_comm_processes, max_no_of_comm_patterns, sync_barrier_mode, &
    & max_mpi_message_size

  USE mo_grid_subset,    ONLY: block_no, index_no
  USE mo_model_domain,    ONLY: t_patch
  USE mo_decomposition_tools, ONLY: t_glb2loc_index_lookup, get_local_index
  USE mo_mpi,             ONLY: p_send, p_irecv, p_wait, p_isend, &
     & p_send, p_real_dp, my_process_is_mpi_seq,   &
     & work_mpi_barrier, &
     & get_my_mpi_work_communicator, get_my_mpi_work_comm_size, &
     & get_my_mpi_work_id
  USE mo_timer,           ONLY: timer_start, timer_stop, timer_icon_comm_sync, &
    & activate_sync_timers, timer_icon_comm_fillrecv, timer_icon_comm_wait, &
    & timer_icon_comm_ircv, timer_icon_comm_fillsend, timer_icon_comm_fillandsend, &
    & timer_icon_comm_isend, timer_icon_comm_barrier_2, timer_icon_comm_send, &
    & timer_start, timer_stop

  USE mo_master_control,  ONLY: get_my_process_name

  USE mo_impl_constants,     ONLY: HALO_LEVELS_CEILING

!   USE mo_impl_constants,  ONLY: &
!     & min_rlcell, max_rlcell, min_rlcell_int, &
!     & min_rlvert, max_rlvert,                 & ! min_rlvert_int,
!     & min_rledge, max_rledge, min_rledge_int


  IMPLICIT NONE

  PRIVATE

  ! public constants

  PUBLIC :: is_ready, until_sync

  ! public methods
  PUBLIC :: construct_icon_comm_lib       ! the first call to the icon_comm_lib
  PUBLIC :: init_icon_std_comm_patterns  ! initilize comm_patterns for a patch
                                     ! must be called before any communication variable is
                                     ! created for this patch
  PUBLIC :: destruct_icon_comm_lib

  PUBLIC :: new_icon_comm_pattern
  PUBLIC :: inverse_of_icon_comm_pattern

  PUBLIC :: new_icon_comm_variable
  PUBLIC :: delete_icon_comm_variable

  PUBLIC :: icon_comm_var_is_ready  ! The comm_variable can be communicated
  PUBLIC :: icon_comm_sync
  PUBLIC :: icon_comm_sync_all

  PUBLIC :: print_grid_comm_pattern, print_grid_comm_stats
!   PUBLIC :: icon_comm_show

  PUBLIC :: t_mpi_mintype
  PUBLIC :: mpi_reduce_mindistance_pts
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! internal parameters
  INTEGER, PARAMETER ::  checksum_size = 2
  ! TAG parameters
  INTEGER, PARAMETER ::  halo_tag = 1
  INTEGER, PARAMETER ::  internal_tag = 2
  !--------------------------------------------------------------
  ! status / requests
  INTEGER, PARAMETER ::  not_active = 0
  INTEGER, PARAMETER ::  active = 1
  INTEGER, PARAMETER ::  is_ready = 10    ! essentially the same as communicate
  INTEGER, PARAMETER ::  communicate = 11

  !--------------------------------------------------------------

  !--------------------------------------------------------------
  ! grid dimensions
  INTEGER, PARAMETER ::  grid_2D = 2
  INTEGER, PARAMETER ::  grid_3D = 3
  !--------------------------------------------------------------
  ! scope parameters
  INTEGER, PARAMETER ::  global = 0
  INTEGER, PARAMETER ::  until_sync = 1

  !--------------------------------------------------------------
  ! TYPE definitions

  !--------------------------------------------------------------
  !> Holds the buffer pointers to procs to be sent and received
  TYPE t_comm_process_buffer

     INTEGER :: pid             ! the process to communicate to

     INTEGER :: start_index      ! the start point in the whole bufffer
     INTEGER :: end_index        ! the end   point in the whole bufffer
     INTEGER :: current_index    ! the current point in the whole bufffer
     INTEGER :: buffer_size      ! the size of this comm to process bufffer

  END TYPE t_comm_process_buffer
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> Holds the communcation pattern for each of the process to which we communicate
  TYPE t_process_comm_pattern

    INTEGER :: pid           ! the process to which we communicate
    INTEGER :: no_of_points  ! the (horizontal) number of values to be communicated
    INTEGER, ALLOCATABLE :: global_index(:) ! the global index of the variable
    INTEGER, ALLOCATABLE :: block_no(:)  ! the local block of the variable
    INTEGER, ALLOCATABLE :: index_no(:)  ! the local index of the variable

    INTEGER :: buffer_index ! the index of the send/recv buffer handler for this process

  END TYPE t_process_comm_pattern
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> Holds the communication patterns to all procs to be communicated.
  !  One structure represents one overall pattern for cells, edges, etc.
  TYPE t_grid_comm_pattern

    INTEGER :: id   ! the id of this, it is the index in the array it lives

    ! the send communication pattern for each process to sent to
    INTEGER :: no_of_send_procs
    TYPE(t_process_comm_pattern), POINTER :: send(:)
    ! the recv communication pattern for each process to receive from
    INTEGER :: no_of_recv_procs
    TYPE(t_process_comm_pattern), POINTER :: recv(:)

    CHARACTER(len=32) :: name

    INTEGER :: status

  END TYPE t_grid_comm_pattern
  !--------------------------------------------------------------

  !--------------------------------------------------------------
  !> Holds the variables to be communicated
  TYPE t_comm_variable_real
    INTEGER :: status                          ! not_active, active, ...
    INTEGER :: request                         ! not_active, ...
    INTEGER :: comm_status                     ! not_active, ...
    INTEGER :: comm_buffer                     ! for future use
    INTEGER :: scope                           ! until_sync = temporary, will be deleted after sync
                                               ! global = will be kept until deleted

    INTEGER :: comm_pattern_index                   ! location where variable lives:
                                                 ! edges, vertices, cells
    INTEGER :: grid_dim                        ! 1D, 2D 3D
    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern    ! the communication pattern

    INTEGER :: no_of_variables                 ! how many variables stored in the array
                                                 ! if =0 then only one

!     INTEGER :: dim_1                         ! should be nproma
    INTEGER :: vertical_layers                 ! if 3D: nlevels
!     INTEGER :: dim_3                         !
    INTEGER :: dim_4                           ! =no_of_variables

    REAL(wp), POINTER :: recv_values_2d(:,:)     ! nproma, nblocks
!    REAL(wp), DIMENSION(:,:,:), POINTER :: recv_values_3d   ! if 3D: nproma, vertical layers, nblocks
    REAL(wp), POINTER :: recv_values_3d(:,:,:)   ! if 3D: nproma, vertical layers, nblocks
    REAL(wp), POINTER :: recv_values_4d(:,:,:,:) ! nproma, vertical layers, nblocks

    REAL(wp), POINTER :: send_values_2d(:,:)     ! nproma, nblocks
!    REAL(wp), DIMENSION(:,:,:), POINTER :: send_values_3d   ! if 3D: nproma, vertical layers, nblocks
    REAL(wp), POINTER :: send_values_3d(:,:,:)   ! if 3D: nproma, vertical layers, nblocks
    REAL(wp), POINTER :: send_values_4d(:,:,:,:) ! nproma, vertical layers, nblocks

    CHARACTER(len=32) :: name

  END TYPE t_comm_variable_real

  !--------------------------------------------------------------

  !> type definition for a user-defined reduction operation, see
  !  subroutine "mpi_reduce_mindistance_pts"
  TYPE t_mpi_mintype
    SEQUENCE
    REAL(wp) :: rdist
    INTEGER  :: glb_index
    INTEGER  :: owner
  END TYPE t_mpi_mintype


  !> The array of the grid objects.
  TYPE(t_comm_variable_real), ALLOCATABLE, TARGET :: comm_variable(:)

  TYPE(t_comm_process_buffer), ALLOCATABLE, TARGET :: send_procs_buffer(:)

  TYPE(t_comm_process_buffer), ALLOCATABLE, TARGET :: recv_procs_buffer(:)

  TYPE(t_grid_comm_pattern), ALLOCATABLE, TARGET :: grid_comm_pattern_list(:)

  ! At the moment we only have one set of communication buffers
  ! This means that only one communation bulk can be executed each time
  REAL(wp), ALLOCATABLE :: send_buffer(:)
  REAL(wp), ALLOCATABLE :: recv_buffer(:)
  INTEGER  :: buffer_comm_status

  !--------------------------------------------------------------
  INTEGER ::  max_comm_patterns = 0
  INTEGER ::  allocated_comm_patterns = 0

  !> The number of actual active comm_variables
  INTEGER :: active_comm_variables
  !> The max id of the active comm_variables
  INTEGER :: max_active_comm_variables

  !> The number of actual active send_buffers
  INTEGER :: active_send_buffers

  !> The number of actual active recv_buffers
  INTEGER :: active_recv_buffers


  LOGICAL :: comm_lib_is_initialized

  INTEGER :: log_file_id = 0

  INTEGER :: max_send_buffer_size, max_recv_buffer_size

  !-------------------------------------------------------------------------
  INTEGER :: my_work_communicator
  INTEGER :: my_work_comm_size
  INTEGER :: my_mpi_work_id
  LOGICAL :: this_is_mpi_sequential

  !-------------------------------------------------------------------------
  INTERFACE new_icon_comm_variable
    MODULE PROCEDURE new_comm_var_r2d
    MODULE PROCEDURE new_comm_var_r2d_recv_send
    MODULE PROCEDURE new_comm_var_r3d
    MODULE PROCEDURE new_comm_var_r3d_recv_send
!     MODULE PROCEDURE new_comm_variable_r3d_target
    MODULE PROCEDURE new_comm_var_r4d
    MODULE PROCEDURE new_comm_var_r4d_recv_send
  END INTERFACE
  !-------------------------------------------------------------------------
  INTERFACE icon_comm_sync
    MODULE PROCEDURE icon_comm_sync_2D_1
    MODULE PROCEDURE icon_comm_sync_3D_1
    MODULE PROCEDURE icon_comm_sync_3D_2
    MODULE PROCEDURE icon_comm_sync_3D_3
    MODULE PROCEDURE icon_comm_sync_3D_4
  END INTERFACE

  !-------------------------------------------------------------------------
CONTAINS

  !-----------------------------------------------------------------------
  !> NOTE: Not completed yet!
  SUBROUTINE destruct_icon_comm_lib()
    CHARACTER(*), PARAMETER :: method_name = "destruct_icon_comm_lib"

    IF (this_is_mpi_sequential) RETURN

    DEALLOCATE(send_buffer, recv_buffer)
    IF (log_file_id > 0) CLOSE(log_file_id)

  END SUBROUTINE destruct_icon_comm_lib
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE construct_icon_comm_lib()

    INTEGER :: i,k, return_status
    LOGICAL :: unit_is_occupied

    CHARACTER(*), PARAMETER :: method_name = "construct_icon_comm_lib"

    this_is_mpi_sequential = my_process_is_mpi_seq()
!     write(0,*) method_name, " this_is_mpi_sequential=", this_is_mpi_sequential

    IF(this_is_mpi_sequential) RETURN


    IF ( comm_lib_is_initialized ) THEN
      CALL message(method_name, 'cannot be called more than once')
      RETURN
    END IF

    ALLOCATE(comm_variable(max_no_of_comm_variables), &
      & stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE comm_variable failed')

    ALLOCATE(send_procs_buffer(max_no_of_comm_processes), &
      & recv_procs_buffer(max_no_of_comm_processes), stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE send-recv_procs_buffer failed')

    allocated_comm_patterns = max_no_of_comm_patterns
    ALLOCATE(grid_comm_pattern_list(allocated_comm_patterns), &
      & stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE grid_comm_pattern_list failed')

    max_send_buffer_size = max_send_recv_buffer_size
    max_recv_buffer_size = max_send_recv_buffer_size
    ALLOCATE(send_buffer(max_send_buffer_size), &
      & recv_buffer(max_recv_buffer_size),stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE send,recv buffers failed')

    active_comm_variables = 0
    max_active_comm_variables = 0
    active_send_buffers = 0
    active_recv_buffers = 0
    comm_lib_is_initialized = .TRUE.
    buffer_comm_status = not_active
    max_comm_patterns = 0

    DO i=1,max_no_of_comm_variables
      CALL clear_comm_variable(i)
    ENDDO

    DO k=1,allocated_comm_patterns
      grid_comm_pattern_list(k)%status = not_active
    ENDDO

    ! initialize the communicators
    my_work_communicator = get_my_mpi_work_communicator()
    my_work_comm_size    = get_my_mpi_work_comm_size()
    my_mpi_work_id       = get_my_mpi_work_id()

    IF (icon_comm_debug) THEN
      DO log_file_id = 500, 5000
        INQUIRE (UNIT=log_file_id, OPENED=unit_is_occupied)
        IF ( .NOT. unit_is_occupied ) EXIT
      ENDDO
      IF (unit_is_occupied) &
        CALL finish(method_name, "Cannot find avaliable file unit")
      WRITE(message_text,'(a,a,a,i4.4)') 'log.', TRIM(get_my_process_name()), &
        & ".icon_comm.", my_mpi_work_id
      OPEN (log_file_id, FILE=TRIM(message_text))
    ENDIF

    RETURN

  END SUBROUTINE construct_icon_comm_lib
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE resize_send_buffer(new_size)
    INTEGER, INTENT(in) :: new_size

    INTEGER :: return_status
    CHARACTER(*), PARAMETER :: method_name = "resize_send_buffer"

    IF (new_size <= max_send_buffer_size) THEN
      CALL warning(method_name, "new_size <= max_send_buffer_size")
      RETURN
    ENDIF

    max_send_buffer_size = new_size + 128
    DEALLOCATE(send_buffer, stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'DEALLOCATE send_buffer failed')
    ALLOCATE(send_buffer(max_send_buffer_size),stat=return_status)
    IF (return_status > 0) &
     CALL finish (method_name, 'ALLOCATE send_buffers failed')

  END SUBROUTINE resize_send_buffer
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE resize_recv_buffer(new_size)
    INTEGER, INTENT(in) :: new_size

    INTEGER :: return_status
    CHARACTER(*), PARAMETER :: method_name = "resize_recv_buffer"

    IF (new_size <= max_recv_buffer_size) THEN
      CALL warning(method_name, "new_size <= max_recv_buffer_size")
      RETURN
    ENDIF

    max_recv_buffer_size = new_size + 128
    DEALLOCATE(recv_buffer, stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'DEALLOCATE recv_buffer failed')
    ALLOCATE(recv_buffer(max_recv_buffer_size),stat=return_status)
    IF (return_status > 0) &
     CALL finish (method_name, 'ALLOCATE recv_buffers failed')

  END SUBROUTINE resize_recv_buffer
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE init_icon_std_comm_patterns(patch)
    TYPE(t_patch), INTENT(inout) :: patch

    CALL init_icon_halo_comm_patterns(patch)

  END SUBROUTINE init_icon_std_comm_patterns
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE init_icon_halo_comm_patterns(patch)
    TYPE(t_patch), INTENT(inout) :: patch

    CHARACTER(*), PARAMETER :: method_name = "init_icon_std_comm_patterns"

    IF(this_is_mpi_sequential) RETURN
!     ! set id of grid_comm_pattern_list to identity
!     DO i = 1, max_comm_patterns
!     ENDDO


    ! halo cells comm_pattern
!     CALL work_mpi_barrier()
!     write(0,*) my_mpi_work_id, method_name, "setup_grid_comm_pattern cells_not_in_domain..."
    patch%sync_cells_not_in_domain = new_icon_comm_pattern( &
      patch%n_patch_cells,   patch%cells%decomp_info%owner_local, &
      patch%cells%decomp_info%glb_index, &
      patch%cells%decomp_info%glb2loc_index, name="cells_not_in_domain" )
    patch%sync_cells_not_owned = patch%sync_cells_not_in_domain

    patch%sync_cells_one_edge_in_domain = new_icon_comm_pattern( &
      patch%n_patch_cells,   patch%cells%decomp_info%owner_local, &
      patch%cells%decomp_info%glb_index, &
      patch%cells%decomp_info%glb2loc_index, &
      halo_level=patch%cells%decomp_info%halo_level, level_start=1, &
      level_end=1, name="cells_one_edge_in_domain" )

    ! halo edges comm_pattern
!     CALL work_mpi_barrier()
!     write(0,*) my_mpi_work_id, method_name, "setup_grid_comm_pattern edges_not_owned..."
    patch%sync_edges_not_owned = new_icon_comm_pattern(   &
      patch%n_patch_edges,   patch%edges%decomp_info%owner_local, &
      patch%edges%decomp_info%glb_index, &
      patch%edges%decomp_info%glb2loc_index, name="edges_not_owned")

    patch%sync_edges_not_in_domain = new_icon_comm_pattern( &
      patch%n_patch_edges,   patch%edges%decomp_info%owner_local, &
      patch%edges%decomp_info%glb_index, &
      patch%edges%decomp_info%glb2loc_index, &
      halo_level=patch%edges%decomp_info%halo_level, level_start=2, &
      level_end=HALO_LEVELS_CEILING, name="edges_not_in_domain")

    ! halo verts comm_pattern
!     CALL work_mpi_barrier()
!     write(0,*) my_mpi_work_id, method_name, "setup_grid_comm_pattern verts_not_owned..."
    patch%sync_verts_not_owned = new_icon_comm_pattern(   &
      patch%n_patch_verts,   patch%verts%decomp_info%owner_local, &
      patch%verts%decomp_info%glb_index, &
      patch%verts%decomp_info%glb2loc_index, name="verts_not_owned" )

    patch%sync_verts_not_in_domain = new_icon_comm_pattern(  &
      patch%n_patch_verts,   patch%verts%decomp_info%owner_local, &
      patch%verts%decomp_info%glb_index, &
      patch%verts%decomp_info%glb2loc_index, &
      halo_level=patch%verts%decomp_info%halo_level, level_start=2, &
      level_end=HALO_LEVELS_CEILING, name="verts_not_in_domain" )

    CALL print_grid_comm_stats(patch%sync_cells_not_in_domain)
    CALL print_grid_comm_stats(patch%sync_cells_one_edge_in_domain)
    CALL print_grid_comm_stats(patch%sync_edges_not_owned)
    CALL print_grid_comm_stats(patch%sync_edges_not_in_domain)
    CALL print_grid_comm_stats(patch%sync_verts_not_owned)
    CALL print_grid_comm_stats(patch%sync_verts_not_in_domain)
    IF ( icon_comm_debug) THEN
      CALL print_grid_comm_pattern(patch%sync_cells_not_in_domain)
      CALL print_grid_comm_pattern(patch%sync_cells_one_edge_in_domain)
      CALL print_grid_comm_pattern(patch%sync_edges_not_owned)
      CALL print_grid_comm_pattern(patch%sync_edges_not_in_domain)
      CALL print_grid_comm_pattern(patch%sync_verts_not_owned)
      CALL print_grid_comm_pattern(patch%sync_verts_not_in_domain)
    ENDIF

 !   CALL finish("init_icon_std_comm_patterns","ends")

  END SUBROUTINE init_icon_halo_comm_patterns
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  INTEGER FUNCTION new_icon_comm_pattern(total_no_of_points, &
    & receive_from_owner, my_global_index, send_glb2loc_index, &
    & allow_send_to_myself, halo_level, level_start, level_end, name)

    INTEGER, INTENT(in) :: total_no_of_points
    INTEGER, INTENT(in) :: receive_from_owner(:), my_global_index(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: send_glb2loc_index
                                          ! global to local index lookup
                                          ! information of the SENDER array
    LOGICAL, INTENT(in), OPTIONAL :: allow_send_to_myself
    INTEGER, INTENT(in), OPTIONAL :: halo_level(:,:), level_start, level_end
    CHARACTER(*), INTENT(in) :: name

    max_comm_patterns = max_comm_patterns + 1
    IF (max_comm_patterns > allocated_comm_patterns) THEN
      CALL finish("new_icon_comm_pattern", &
        &         "max_comm_patterns > allocated_comm_patterns")
    ENDIF

    grid_comm_pattern_list(max_comm_patterns)%id = max_comm_patterns
    new_icon_comm_pattern = max_comm_patterns

    CALL setup_grid_comm_pattern(grid_comm_pattern_list(max_comm_patterns), &
      & total_no_of_points, receive_from_owner, my_global_index,            &
      & send_glb2loc_index, allow_send_to_myself, halo_level, level_start,    &
      & level_end, name)


  END FUNCTION new_icon_comm_pattern
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE setup_grid_comm_pattern(grid_comm_pattern, total_no_of_points, &
    & receive_from_owner, my_global_index, send_glb2loc_index,  allow_send_to_myself,   &
    & halo_level, level_start, level_end, name)

    TYPE(t_grid_comm_pattern), INTENT(inout) :: grid_comm_pattern
    INTEGER, INTENT(in) :: total_no_of_points
    INTEGER, INTENT(in) :: receive_from_owner(:), my_global_index(:)
    TYPE(t_glb2loc_index_lookup), INTENT(IN) :: send_glb2loc_index
                                          ! global to local index lookup
                                          ! information of the SENDER array
    LOGICAL, INTENT(in), OPTIONAL :: allow_send_to_myself
    INTEGER, INTENT(in), OPTIONAL :: halo_level(:,:), level_start, level_end
    CHARACTER(*), INTENT(in) :: name

    TYPE(t_process_comm_pattern), POINTER :: p_comm_pattern

    INTEGER :: comm_points(max_no_of_comm_processes)
    INTEGER :: procs_id(max_no_of_comm_processes)
    INTEGER :: buffer_id(max_no_of_comm_processes)
    INTEGER :: comm_of_buffer_id(max_no_of_comm_processes)
    INTEGER :: filtered_receive_from_owner(total_no_of_points)

    INTEGER, ALLOCATABLE :: recv_requests(:), total_requests(:)
    INTEGER, ALLOCATABLE :: recv_global_indexes(:,:), send_global_indexes(:,:)

    INTEGER :: this_mpi_work_id
    INTEGER :: i, point_idx, bfid, no_comm_procs, no_of_points
    INTEGER :: owner_id, max_comm_points, max_buffer_size, no_of_recv_requests
    INTEGER :: local_idx, return_status

    CHARACTER(*), PARAMETER :: method_name = "setup_grid_comm_pattern"


  END SUBROUTINE setup_grid_comm_pattern
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  INTEGER FUNCTION inverse_of_icon_comm_pattern(in_comm_pattern_id, name)

    INTEGER, INTENT(in) :: in_comm_pattern_id
    CHARACTER(*), INTENT(in) :: name

    TYPE(t_grid_comm_pattern), POINTER :: initial_comm_pattern, inverse_comm_pattern
    TYPE(t_process_comm_pattern), POINTER :: inv_proc_comm_pattern
    INTEGER :: no_of_recv_procs, no_of_send_procs, return_status, i, j, no_of_points

    CHARACTER(*), PARAMETER :: method_name = "inverse_of_icon_comm_pattern"

    ! get next avail comm_patterns space
    max_comm_patterns = max_comm_patterns + 1
    IF (max_comm_patterns > allocated_comm_patterns) THEN
      CALL finish("new_icon_comm_pattern","max_comm_patterns > allocated_comm_patterns")
    ENDIF

    initial_comm_pattern => grid_comm_pattern_list(in_comm_pattern_id)
    inverse_comm_pattern => grid_comm_pattern_list(max_comm_patterns)

    inverse_comm_pattern%id = max_comm_patterns

    no_of_recv_procs = initial_comm_pattern%no_of_send_procs
    no_of_send_procs = initial_comm_pattern%no_of_recv_procs
    inverse_comm_pattern%no_of_recv_procs = no_of_recv_procs
    inverse_comm_pattern%no_of_send_procs = no_of_send_procs

    ALLOCATE(inverse_comm_pattern%recv(no_of_recv_procs),stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE(inverse_comm_pattern%recv()')
    ALLOCATE(inverse_comm_pattern%send(no_of_send_procs),stat=return_status)
    IF (return_status > 0) &
      CALL finish (method_name, 'ALLOCATE(inverse_comm_pattern%send()')

    ! fill recv patterns from the initial sent patterns
    DO i=1, no_of_recv_procs
      inv_proc_comm_pattern => inverse_comm_pattern%recv(i)
      no_of_points = initial_comm_pattern%send(i)%no_of_points

      inv_proc_comm_pattern%pid = initial_comm_pattern%send(i)%pid
      inv_proc_comm_pattern%no_of_points = no_of_points
      inv_proc_comm_pattern%buffer_index = &
        & get_recvbuffer_id_of_pid(inv_proc_comm_pattern%pid)

      ALLOCATE(inv_proc_comm_pattern%block_no(no_of_points),     &
        &      inv_proc_comm_pattern%index_no(no_of_points),     &
        &      inv_proc_comm_pattern%global_index(no_of_points), &
        &      stat=return_status)
      IF (return_status > 0) &
        CALL finish (method_name, 'ALLOCATE(process_comm_pattern arrays)')

      DO j=1, no_of_points
        inv_proc_comm_pattern%block_no(j)     = initial_comm_pattern%send(i)%block_no(j)
        inv_proc_comm_pattern%index_no(j)     = initial_comm_pattern%send(i)%index_no(j)
        inv_proc_comm_pattern%global_index(j) = initial_comm_pattern%send(i)%global_index(j)
      ENDDO

    ENDDO ! i=1, no_of_recv_procs

    ! fill send patterns from the initial receive patterns
    DO i=1, no_of_send_procs
      inv_proc_comm_pattern => inverse_comm_pattern%send(i)
      no_of_points = initial_comm_pattern%recv(i)%no_of_points

      inv_proc_comm_pattern%pid = initial_comm_pattern%recv(i)%pid
      inv_proc_comm_pattern%no_of_points = no_of_points
      inv_proc_comm_pattern%buffer_index = &
        & get_sendbuffer_id_of_pid(inv_proc_comm_pattern%pid)

      ALLOCATE(inv_proc_comm_pattern%block_no(no_of_points),     &
        &      inv_proc_comm_pattern%index_no(no_of_points),     &
        &      inv_proc_comm_pattern%global_index(no_of_points), &
        &      stat=return_status)
      IF (return_status > 0) &
        CALL finish (method_name, 'ALLOCATE(process_comm_pattern arrays)')

      DO j=1, no_of_points
        inv_proc_comm_pattern%block_no(j)     = initial_comm_pattern%recv(i)%block_no(j)
        inv_proc_comm_pattern%index_no(j)     = initial_comm_pattern%recv(i)%index_no(j)
        inv_proc_comm_pattern%global_index(j) = initial_comm_pattern%recv(i)%global_index(j)
      ENDDO

    ENDDO ! i=1, no_of_send_procs
   !--------------------------------------------------------------

    inverse_comm_pattern%status = active
    inverse_comm_pattern%name   = TRIM(name)
    inverse_of_icon_comm_pattern = max_comm_patterns

  END FUNCTION inverse_of_icon_comm_pattern
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE print_grid_comm_stats(comm_pattern_id)

    INTEGER, INTENT(in) :: comm_pattern_id

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    INTEGER :: i, min_points, max_points, tot_points

    IF ( log_file_id <= 0 ) RETURN

    grid_comm_pattern => grid_comm_pattern_list(comm_pattern_id)
    write(log_file_id,*) " === Communication stats for ", TRIM(grid_comm_pattern%name), &
      & " id=", grid_comm_pattern%id
    write(log_file_id,*) " --- recv: no_of_procs= ", grid_comm_pattern%no_of_recv_procs
    min_points = 99999999
    max_points = 0
    tot_points = 0
    DO i=1,grid_comm_pattern%no_of_recv_procs
      min_points = MIN(min_points, grid_comm_pattern%recv(i)%no_of_points)
      max_points = MAX(max_points, grid_comm_pattern%recv(i)%no_of_points)
      tot_points = tot_points + grid_comm_pattern%recv(i)%no_of_points
    ENDDO
    write(log_file_id,*) " - min halos=",min_points, " max halos=",max_points, &
      & " total halos=", tot_points

    write(log_file_id,*) " --- send: no_of_procs= ", grid_comm_pattern%no_of_send_procs
    min_points = 99999999
    max_points = 0
    tot_points = 0
    DO i=1,grid_comm_pattern%no_of_send_procs
      min_points = MIN(min_points, grid_comm_pattern%send(i)%no_of_points)
      max_points = MAX(max_points, grid_comm_pattern%send(i)%no_of_points)
      tot_points = tot_points + grid_comm_pattern%send(i)%no_of_points
    ENDDO
    write(log_file_id,*) " - min halos=",min_points, " max halos=",max_points, &
      & " total halos=", tot_points


    write(log_file_id,*) " === End Communication stats for ", TRIM(grid_comm_pattern%name)

  END SUBROUTINE print_grid_comm_stats
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  SUBROUTINE print_grid_comm_pattern(comm_pattern_id)
    INTEGER, INTENT(in) :: comm_pattern_id

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    INTEGER :: i

    IF ( log_file_id <= 0 ) RETURN

    grid_comm_pattern => grid_comm_pattern_list(comm_pattern_id)

    write(log_file_id,*) " === Communication pattern info for ", TRIM(grid_comm_pattern%name), &
      & " id=", grid_comm_pattern%id
    write(log_file_id,*) " - no_of_recv_procs= ", grid_comm_pattern%no_of_recv_procs
    write(log_file_id,*) " - no_of_send_procs= ", grid_comm_pattern%no_of_send_procs

    DO i=1,grid_comm_pattern%no_of_recv_procs
      write(log_file_id,*) " - recv from ", grid_comm_pattern%recv(i)%pid, " size=", &
        & grid_comm_pattern%recv(i)%no_of_points
      write(log_file_id,*) " - global indexes=", grid_comm_pattern%recv(i)%global_index(:)
    ENDDO

    DO i=1,grid_comm_pattern%no_of_send_procs
      write(log_file_id,*) " - send to ", grid_comm_pattern%send(i)%pid, " size=", &
        & grid_comm_pattern%send(i)%no_of_points
      write(log_file_id,*) " - global indexes=", grid_comm_pattern%send(i)%global_index(:)
    ENDDO

    write(log_file_id,*) " === End Communication pattern info for ", TRIM(grid_comm_pattern%name)

  END SUBROUTINE print_grid_comm_pattern
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  INTEGER FUNCTION get_recvbuffer_id_of_pid( pid )
    INTEGER, INTENT(IN)       :: pid

    INTEGER :: bfid
    CHARACTER(*), PARAMETER :: method_name = "get_recvbuffer_id_of_pid"

    DO bfid = 1, active_recv_buffers
      IF ( recv_procs_buffer(bfid)%pid == pid ) THEN
        get_recvbuffer_id_of_pid = bfid
        RETURN
      ENDIF
    ENDDO

    ! add a new recv_procs_buffer
    IF (active_recv_buffers >= max_no_of_comm_processes ) &
      CALL finish(method_name, "active_recv_buffers >= max_no_of_comm_processes")
    active_recv_buffers = active_recv_buffers + 1
    recv_procs_buffer(active_recv_buffers)%pid = pid
    get_recvbuffer_id_of_pid = active_recv_buffers

  END FUNCTION get_recvbuffer_id_of_pid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  INTEGER FUNCTION get_sendbuffer_id_of_pid( pid )
    INTEGER, INTENT(IN)       :: pid

    INTEGER :: bfid
    CHARACTER(*), PARAMETER :: method_name = "get_sendbuffer_id_of_pid"

    DO bfid = 1, active_send_buffers
      IF ( send_procs_buffer(bfid)%pid == pid ) THEN
        get_sendbuffer_id_of_pid = bfid
        RETURN
      ENDIF
    ENDDO

    ! add a new send_procs_buffer
    IF (active_send_buffers >= max_no_of_comm_processes ) &
      CALL finish(method_name, "active_send_buffers >= max_no_of_comm_processes")
    active_send_buffers = active_send_buffers + 1
    send_procs_buffer(active_send_buffers)%pid = pid
    get_sendbuffer_id_of_pid = active_send_buffers

  END FUNCTION get_sendbuffer_id_of_pid
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_var_r4d(var, comm_pattern_index, vertical_layers, &
    & no_of_variables, status, scope, name )
!     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:,:)
    REAL(wp), TARGET   :: var(:,:,:,:)
    INTEGER, INTENT(IN) :: comm_pattern_index

     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
!     INTEGER, INTENT(IN), OPTIONAL :: var_dim
    INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
    INTEGER, INTENT(IN), OPTIONAL :: status
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    new_comm_var_r4d = new_comm_var_r4d_recv_send( &
      & var, var, comm_pattern_index, vertical_layers, &
      & no_of_variables, status, scope, name )

  END FUNCTION new_comm_var_r4d
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_var_r4d_recv_send(recv_var, send_var, comm_pattern_index, &
    & vertical_layers, no_of_variables, status, scope, name )
!     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:,:)
    REAL(wp), TARGET  :: recv_var(:,:,:,:)
    REAL(wp), TARGET  :: send_var(:,:,:,:)
    INTEGER, INTENT(IN)       :: comm_pattern_index

     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
!     INTEGER, INTENT(IN), OPTIONAL :: var_dim
    INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
    INTEGER, INTENT(IN), OPTIONAL :: status
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    REAL(wp), POINTER  :: recv_var_3d(:,:,:), send_var_3d(:,:,:)

    CHARACTER(*), PARAMETER :: method_name = "new_comm_var_r4d"


    IF  (this_is_mpi_sequential) THEN
      new_comm_var_r4d_recv_send = 0
      RETURN
    ENDIF

    recv_var_3d => recv_var(:,:,:,1)
    send_var_3d => send_var(:,:,:,1)

    IF (PRESENT(vertical_layers)) THEN
      new_comm_var_r4d_recv_send = &
        & new_comm_var_r3d_recv_send(recv_var=recv_var_3d, send_var=send_var_3d, &
        & comm_pattern_index = comm_pattern_index, vertical_layers = vertical_layers)
    ELSE
      new_comm_var_r4d_recv_send = &
        & new_comm_var_r3d_recv_send(recv_var=recv_var_3d, send_var=send_var_3d, &
        & comm_pattern_index = comm_pattern_index, vertical_layers = vertical_layers)
    ENDIF

    comm_variable(new_comm_var_r4d_recv_send)%recv_values_4d => recv_var
    comm_variable(new_comm_var_r4d_recv_send)%send_values_4d => send_var

    NULLIFY(comm_variable(new_comm_var_r4d_recv_send)%recv_values_3d)
    NULLIFY(comm_variable(new_comm_var_r4d_recv_send)%send_values_3d)

    IF (PRESENT(no_of_variables)) THEN
      comm_variable(new_comm_var_r4d_recv_send)%no_of_variables = no_of_variables
    ELSE
      comm_variable(new_comm_var_r4d_recv_send)%no_of_variables = SIZE(recv_var,4)
    ENDIF

    comm_variable(new_comm_var_r4d_recv_send)%dim_4 = &
      comm_variable(new_comm_var_r4d_recv_send)%no_of_variables

    IF (PRESENT(status)) THEN
      IF (status == is_ready) CALL icon_comm_var_is_ready(new_comm_var_r4d_recv_send)
    ENDIF

    IF (PRESENT(scope)) THEN
      comm_variable(new_comm_var_r4d_recv_send)%scope = scope
    ELSE
      comm_variable(new_comm_var_r4d_recv_send)%scope = global
    ENDIF

    IF (PRESENT(name)) THEN
      comm_variable(new_comm_var_r4d_recv_send)%name = TRIM(name)
    ELSE
      comm_variable(new_comm_var_r4d_recv_send)%name = ""
    ENDIF

  END FUNCTION new_comm_var_r4d_recv_send
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_var_r3d(var, comm_pattern_index,  &
    & vertical_layers, status, scope, name)

    REAL(wp), TARGET   :: var(:,:,:)
    INTEGER, INTENT(IN) :: comm_pattern_index

    INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
    INTEGER, INTENT(IN), OPTIONAL :: status
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    new_comm_var_r3d =  new_comm_var_r3d_recv_send( &
      & var, var, comm_pattern_index,  &
      & vertical_layers, status, scope, name)

  END FUNCTION new_comm_var_r3d
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_var_r3d_recv_send(recv_var, send_var, comm_pattern_index,  &
    & vertical_layers, status, scope, name) !, &
   ! & var_dim, no_of_variables, vertical_layers)

!     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:)
    REAL(wp), TARGET   :: recv_var(:,:,:)
    REAL(wp), TARGET   :: send_var(:,:,:)
    INTEGER, INTENT(IN) :: comm_pattern_index

!     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
!     INTEGER, INTENT(IN), OPTIONAL :: var_dim
    INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
    INTEGER, INTENT(IN), OPTIONAL :: status
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    CHARACTER(*), PARAMETER :: method_name = "new_comm_var_r3d"

    IF(this_is_mpi_sequential) THEN
      new_comm_var_r3d_recv_send = 0
      RETURN
    ENDIF

    new_comm_var_r3d_recv_send = get_new_comm_variable()

    comm_variable(new_comm_var_r3d_recv_send)%comm_pattern_index = comm_pattern_index

    ! check if comm_pattern is initialized
    IF ( grid_comm_pattern_list(comm_pattern_index)%status == not_active ) &
      CALL finish(method_name, "grid_comm_pattern status = not_active")

    comm_variable(new_comm_var_r3d_recv_send)%grid_comm_pattern => &
        & grid_comm_pattern_list(comm_pattern_index)

    ! this is for a 3D variable
    comm_variable(new_comm_var_r3d_recv_send)%grid_dim = grid_3D

    ! check the dim_1
!     IF( SIZE(var,1) /= nproma ) THEN
!        CALL finish(method_name, 'SIZE(var,1) /= nproma')
!     ENDIF
!     comm_variable(new_comm_var_r3d)%dim_1 = SIZE(var,1)

    ! check the vertical_layers
    comm_variable(new_comm_var_r3d_recv_send)%vertical_layers =  &
      & SIZE(recv_var,2)
    IF ( PRESENT(vertical_layers) ) THEN
      IF ( vertical_layers <=  &
        & comm_variable(new_comm_var_r3d_recv_send)%vertical_layers) THEN
         comm_variable(new_comm_var_r3d_recv_send)%vertical_layers = vertical_layers
      ELSE
         CALL finish(method_name, "vertical_layers are greater than SIZE(var,2)")
      END IF
    END IF

    comm_variable(new_comm_var_r3d_recv_send)%recv_values_3d => recv_var
    comm_variable(new_comm_var_r3d_recv_send)%send_values_3d => send_var
    comm_variable(new_comm_var_r3d_recv_send)%no_of_variables = 1

    IF (PRESENT(status)) THEN
      IF (status == is_ready) CALL icon_comm_var_is_ready(new_comm_var_r3d_recv_send )
    ENDIF

    IF (PRESENT(scope)) THEN
      comm_variable(new_comm_var_r3d_recv_send)%scope = scope
    ELSE
      comm_variable(new_comm_var_r3d_recv_send)%scope = global
    ENDIF

    IF (PRESENT(name)) THEN
      comm_variable(new_comm_var_r3d_recv_send)%name = TRIM(name)
    ELSE
      comm_variable(new_comm_var_r3d_recv_send)%name = ""
    ENDIF

  END FUNCTION new_comm_var_r3d_recv_send
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_var_r2d(var, comm_pattern_index, &
    & vertical_layers, status, scope, name) !, &
   ! & var_dim, no_of_variables, vertical_layers)

!     REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:)
    REAL(wp), TARGET   :: var(:,:)

    INTEGER, INTENT(IN)       :: comm_pattern_index

!     INTEGER, INTENT(IN), OPTIONAL :: no_of_variables
!     INTEGER, INTENT(IN), OPTIONAL :: var_dim
     INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     TYPE(t_comm_pattern), INTENT(IN), POINTER, OPTIONAL :: comm_pattern
    INTEGER, INTENT(IN), OPTIONAL :: status
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    new_comm_var_r2d = new_comm_var_r2d_recv_send( &
      & var, var, comm_pattern_index, &
      & vertical_layers, status, scope, name)

  END FUNCTION new_comm_var_r2d
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
!   INTEGER FUNCTION new_comm_var_r2d_recv_send_t(target_flag, recv_var, send_var, comm_pattern_index, &
!     & vertical_layers, status, scope, name) !, &
!
!     CHARACTER, INTENT(IN) :: target_flag
!     REAL(wp), TARGET   :: recv_var(:,:)
!     REAL(wp), TARGET   :: send_var(:,:)
!
!     INTEGER, INTENT(IN) :: comm_pattern_index
!
!     INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
!     INTEGER, INTENT(IN), OPTIONAL :: status
!     INTEGER, INTENT(IN), OPTIONAL :: scope
!     CHARACTER(*), INTENT(IN), OPTIONAL :: name
!
!     CALL new_comm_var_r2d_recv_send(recv_var, send_var, comm_pattern_index, &
!     & vertical_layers, status, scope, name)
!
!    END FUNCTION new_comm_var_r2d_recv_send_t
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Creates a new comm_variable and returns its id.
  INTEGER FUNCTION new_comm_var_r2d_recv_send(recv_var, send_var, comm_pattern_index, &
    & vertical_layers, status, scope, name) 

    REAL(wp), TARGET   :: recv_var(:,:)
    REAL(wp), TARGET   :: send_var(:,:)
!     REAL(wp), POINTER   :: recv_var(:,:)
!     REAL(wp), POINTER   :: send_var(:,:)

    INTEGER, INTENT(IN) :: comm_pattern_index

    INTEGER, INTENT(IN), OPTIONAL :: vertical_layers
    INTEGER, INTENT(IN), OPTIONAL :: status
    INTEGER, INTENT(IN), OPTIONAL :: scope
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    CHARACTER(*), PARAMETER :: method_name = "new_comm_var_r2d"

    IF(this_is_mpi_sequential) THEN
      new_comm_var_r2d_recv_send = 0
      RETURN
    ENDIF

    new_comm_var_r2d_recv_send = get_new_comm_variable()

    ! where the variable lives
!     SELECT CASE (comm_pattern_index)
!     CASE ( cells_not_in_domain, cells_one_edge_in_domain, edges_not_owned, verts_not_owned)
!       comm_variable(new_comm_var_r2d)%comm_pattern_index = comm_pattern_index
!     CASE default
!       CALL finish(method_name, "Unrecoginzed comm_pattern_index")
!     END SELECT
    comm_variable(new_comm_var_r2d_recv_send)%comm_pattern_index = comm_pattern_index

    ! check if comm_pattern is initialized
    IF ( grid_comm_pattern_list(comm_pattern_index)%status == not_active ) &
      CALL finish(method_name, "grid_comm_pattern status = not_active")

    comm_variable(new_comm_var_r2d_recv_send)%grid_comm_pattern => &
        & grid_comm_pattern_list(comm_pattern_index)

    ! this is for a 3D variable
    comm_variable(new_comm_var_r2d_recv_send)%grid_dim = grid_2D

    ! check the dim_1
!     IF( SIZE(var,1) /= nproma ) THEN
!        CALL finish(method_name, 'SIZE(var,1) /= nproma')
!     ENDIF
!     comm_variable(new_comm_var_r2d)%dim_1 = nproma

    ! check the vertical_layers
    IF (PRESENT(vertical_layers)) THEN
      IF (vertical_layers /= 1) THEN
        CALL finish(method_name, "Unsupported number of layers.")
      ENDIF
    ENDIF
    comm_variable(new_comm_var_r2d_recv_send)%vertical_layers = 1

    comm_variable(new_comm_var_r2d_recv_send)%recv_values_2d => recv_var
    comm_variable(new_comm_var_r2d_recv_send)%send_values_2d => send_var
    comm_variable(new_comm_var_r2d_recv_send)%no_of_variables = 1

    IF (PRESENT(status)) THEN
      IF (status == is_ready) CALL icon_comm_var_is_ready(new_comm_var_r2d_recv_send )
    ENDIF

    IF (PRESENT(scope)) THEN
      comm_variable(new_comm_var_r2d_recv_send)%scope = scope
    ELSE
      comm_variable(new_comm_var_r2d_recv_send)%scope = global
    ENDIF

    IF (PRESENT(name)) THEN
      comm_variable(new_comm_var_r2d_recv_send)%name = TRIM(name)
    ELSE
      comm_variable(new_comm_var_r2d_recv_send)%name = ""
    ENDIF

  END FUNCTION new_comm_var_r2d_recv_send
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! Returns the id of a new comm_variable
  INTEGER FUNCTION get_new_comm_variable()

    INTEGER :: i

    CHARACTER(*), PARAMETER :: method_name = "get_new_comm_variable"


    IF ( .NOT. comm_lib_is_initialized ) &
      CALL finish(method_name, 'comm_lib is not initialized')

    IF (max_active_comm_variables /= active_comm_variables) THEN
      ! we have a hole (inactive grid) in the list of max_active_grids
      DO i = 1, max_active_comm_variables
        IF ( comm_variable(i)%status == not_active ) THEN
          get_new_comm_variable = i
          EXIT
        ENDIF
      ENDDO
    ELSE
      ! check if we reached the end
      IF (active_comm_variables >= max_no_of_comm_variables) THEN
        CALL finish(method_name, 'exceeded max_no_of_comm_variables');
      ENDIF
      max_active_comm_variables = max_active_comm_variables + 1
      get_new_comm_variable = max_active_comm_variables
    ENDIF

    active_comm_variables = active_comm_variables + 1

    comm_variable(get_new_comm_variable)%status = active

  END FUNCTION get_new_comm_variable
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Deletes the comm_variable
  SUBROUTINE delete_icon_comm_variable(comm_variable_id)
    INTEGER, INTENT(in) :: comm_variable_id

    CHARACTER(*), PARAMETER :: method_name = "delete_icon_comm_variable"

!     CALL check_active_comm_variable(comm_variable_id)
    IF(this_is_mpi_sequential) RETURN

    comm_variable(comm_variable_id)%status = not_active
    active_comm_variables = active_comm_variables - 1
    IF ( comm_variable_id == max_active_comm_variables ) &
      max_active_comm_variables = max_active_comm_variables - 1

    CALL clear_comm_variable(comm_variable_id)

  END SUBROUTINE delete_icon_comm_variable
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE clear_comm_variable(id)
    INTEGER, INTENT(in) :: id

    comm_variable(id)%status = not_active
    comm_variable(id)%request = not_active
    comm_variable(id)%comm_status = not_active
    comm_variable(id)%comm_buffer = 0

    comm_variable(id)%no_of_variables = 0
!     comm_variable(id)%dim_1 = 0
    comm_variable(id)%vertical_layers = 0
!     comm_variable(id)%dim_3 = 0
    comm_variable(id)%dim_4 = 0

    NULLIFY(comm_variable(id)%recv_values_2d)
    NULLIFY(comm_variable(id)%send_values_2d)
    NULLIFY(comm_variable(id)%recv_values_3d)
    NULLIFY(comm_variable(id)%send_values_3d)
    NULLIFY(comm_variable(id)%recv_values_4d)
    NULLIFY(comm_variable(id)%send_values_4d)

  END SUBROUTINE clear_comm_variable
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  !! The comm_variable can be communicated
  SUBROUTINE icon_comm_var_is_ready(comm_variable_id)
    INTEGER, INTENT(in) :: comm_variable_id
    TYPE(t_comm_variable_real),pointer :: var

     IF(this_is_mpi_sequential) RETURN

!     CALL check_active_comm_variable(comm_variable_id)
!$OMP SINGLE
      var => comm_variable(comm_variable_id)
      var%request = communicate
!$OMP END SINGLE

  END SUBROUTINE icon_comm_var_is_ready
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_2D_1(var, comm_pattern_index, name)
    INTEGER, INTENT(IN)       :: comm_pattern_index
!    REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:)
    REAL(wp), TARGET   :: var(:,:)
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    INTEGER :: comm_var

    IF(this_is_mpi_sequential) RETURN

    comm_var = new_icon_comm_variable(var,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name )
    CALL icon_comm_sync_all

  END SUBROUTINE icon_comm_sync_2D_1
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_3D_1(var,  comm_pattern_index, name)
    INTEGER, INTENT(IN)       :: comm_pattern_index
!    REAL(wp), POINTER, INTENT(INOUT)   :: var(:,:,:)
    REAL(wp), TARGET   :: var(:,:,:)
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    INTEGER :: comm_var

    IF(this_is_mpi_sequential) RETURN

    comm_var = new_icon_comm_variable(var,  comm_pattern_index,  &
      & status=is_ready, scope=until_sync, name=name )
    CALL icon_comm_sync_all

  END SUBROUTINE icon_comm_sync_3D_1
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_3D_2(var1,  var2, comm_pattern_index, name)
    INTEGER, INTENT(IN)       :: comm_pattern_index
!     REAL(wp), POINTER, INTENT(INOUT)   :: var1(:,:,:)
!     REAL(wp), POINTER, INTENT(INOUT)   :: var2(:,:,:)
    REAL(wp), TARGET   :: var1(:,:,:)
    REAL(wp), TARGET   :: var2(:,:,:)
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    INTEGER :: comm_var_1, comm_var_2

    IF(this_is_mpi_sequential) RETURN

    comm_var_1 = new_icon_comm_variable(var1,  comm_pattern_index,  &
      & status=is_ready, scope=until_sync, name=name)
    comm_var_2 = new_icon_comm_variable(var2,  comm_pattern_index,  &
      & status=is_ready, scope=until_sync, name=name)
    CALL icon_comm_sync_all

  END SUBROUTINE icon_comm_sync_3D_2
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_3D_3(var1,  var2, var3, comm_pattern_index, name)
    INTEGER, INTENT(IN)       :: comm_pattern_index
!     REAL(wp), POINTER, INTENT(INOUT)   :: var1(:,:,:)
!     REAL(wp), POINTER, INTENT(INOUT)   :: var2(:,:,:)
    REAL(wp), TARGET   :: var1(:,:,:)
    REAL(wp), TARGET   :: var2(:,:,:)
    REAL(wp), TARGET   :: var3(:,:,:)
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    INTEGER :: comm_var_1, comm_var_2, comm_var_3

    IF(this_is_mpi_sequential) RETURN

    comm_var_1 = new_icon_comm_variable(var1,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name)
    comm_var_2 = new_icon_comm_variable(var2,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name)
    comm_var_3 = new_icon_comm_variable(var3,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name)
    CALL icon_comm_sync_all

  END SUBROUTINE icon_comm_sync_3D_3
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE icon_comm_sync_3D_4(var1,  var2, var3, var4, comm_pattern_index, name)
    INTEGER, INTENT(IN)       :: comm_pattern_index
!     REAL(wp), POINTER, INTENT(INOUT)   :: var1(:,:,:)
!     REAL(wp), POINTER, INTENT(INOUT)   :: var2(:,:,:)
    REAL(wp), TARGET   :: var1(:,:,:)
    REAL(wp), TARGET   :: var2(:,:,:)
    REAL(wp), TARGET   :: var3(:,:,:)
    REAL(wp), TARGET   :: var4(:,:,:)
    CHARACTER(*), INTENT(IN), OPTIONAL :: name

    INTEGER :: comm_var_1, comm_var_2, comm_var_3, comm_var_4

    IF(this_is_mpi_sequential) RETURN

    comm_var_1 = new_icon_comm_variable(var1,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name)
    comm_var_2 = new_icon_comm_variable(var2,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name)
    comm_var_3 = new_icon_comm_variable(var3,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name)
    comm_var_4 = new_icon_comm_variable(var4,  comm_pattern_index, &
      & status=is_ready, scope=until_sync, name=name)
    CALL icon_comm_sync_all

  END SUBROUTINE icon_comm_sync_3D_4
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE icon_comm_sync_all()

    LOGICAL :: exist_communication_var
    INTEGER :: comm_var

    CHARACTER(*), PARAMETER :: method_name = "icon_comm_sync_all"



    IF(this_is_mpi_sequential) RETURN

    ! check if we have any variables to communicate
    exist_communication_var = .false.
    DO comm_var = 1, max_active_comm_variables
      IF ( comm_variable(comm_var)%request == communicate ) THEN
        exist_communication_var = .true.
        EXIT
      ENDIF
    ENDDO

    IF (.NOT. exist_communication_var) RETURN

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_sync)

    IF (buffer_comm_status == not_active) THEN
      ! no communication steps have been taken
      ! go through the whole process
      SELECT CASE(icon_comm_method)

      CASE(1,101)
        CALL nonblockrecv_all_data()
        CALL fill_send_buffers()
        CALL nonblocksent_all_data()
        ! Wait for all outstanding requests to finish
        IF (activate_sync_timers) CALL timer_start(timer_icon_comm_wait)
        CALL p_wait
        IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_wait)
        CALL fill_vars_from_recv_buffers()

      CASE(2,102)
        CALL nonblockrecv_all_data()
        CALL fill_and_nonblocksend_buffers()
        ! Wait for all outstanding requests to finish
        IF (activate_sync_timers) CALL timer_start(timer_icon_comm_wait)
        CALL p_wait
        IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_wait)
        CALL fill_vars_from_recv_buffers()

      CASE(3,103)
        CALL nonblockrecv_all_data()
        CALL fill_send_buffers()
        CALL blocksent_all_data()
        ! Wait for all outstanding requests to finish
        IF (activate_sync_timers) CALL timer_start(timer_icon_comm_wait)
        CALL p_wait
        IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_wait)
        CALL fill_vars_from_recv_buffers()

      CASE DEFAULT
        CALL finish( method_name,'unknown icon_comm_method.')
      END SELECT

      CALL clear_comm()

    ELSE

      CALL finish(method_name, 'buffer_comm_status /= not_active');

    ENDIF

    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_sync)

    IF (sync_barrier_mode == 2) THEN
      CALL timer_start(timer_icon_comm_barrier_2)
      CALL work_mpi_barrier()
      CALL timer_stop(timer_icon_comm_barrier_2)
   ENDIF

  END SUBROUTINE icon_comm_sync_all
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE clear_comm()

    INTEGER :: comm_var! ,  bfid

    buffer_comm_status = not_active

    IF(this_is_mpi_sequential) RETURN

    ! clear buffer sizes
!     DO bfid = 1, active_send_buffers
!       send_procs_buffer(bfid)%buffer_size = 0
!     ENDDO

    ! go through all variables and clean comm requests
    DO comm_var = 1, max_active_comm_variables
      IF (comm_variable(comm_var)%scope == until_sync) THEN
         CALL delete_icon_comm_variable(comm_var)
      ELSE
        comm_variable(comm_var)%request     = not_active
        comm_variable(comm_var)%comm_status = not_active
      ENDIF
    ENDDO

    CALL adjust_max_active_variables()

  END SUBROUTINE clear_comm
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE adjust_max_active_variables()
    INTEGER :: comm_var, new_max_active_variables

    new_max_active_variables=0
    DO comm_var = 1, max_active_comm_variables
      IF ( comm_variable(comm_var)%status /= not_active ) &
        new_max_active_variables = comm_var
    ENDDO
    max_active_comm_variables = new_max_active_variables

  END SUBROUTINE adjust_max_active_variables
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  SUBROUTINE fill_send_buffers()

    INTEGER :: comm_var, var_no
!     CHARACTER(*), PARAMETER :: method_name = "compute_send_buffer_sizes"

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_fillsend)

    CALL compute_send_buffer_sizes()

    DO comm_var = 1, max_active_comm_variables
      !
      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE

        DO var_no = 1, comm_variable(comm_var)%no_of_variables
          CALL fill_send_buffers_var(comm_var, var_no)
        END DO

    END DO
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_fillsend)

  END SUBROUTINE fill_send_buffers
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE fill_and_nonblocksend_buffers()

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    INTEGER :: comm_var, var_no, bfid, np, buffer_start, buffer_size,&
      & message_size, message_seq_id

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_fillandsend)

    CALL compute_send_buffer_sizes()

!ICON_OMP_PARALLEL IF(icon_comm_openmp)
!ICON_OMP_DO PRIVATE(bfid, comm_var, grid_comm_pattern, np, var_no)
    DO bfid = 1, active_send_buffers
      IF ( send_procs_buffer(bfid)%buffer_size == 0) CYCLE

      DO comm_var = 1, max_active_comm_variables
        IF ( comm_variable(comm_var)%request /= communicate ) CYCLE
        grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
!         vertical_layers = comm_variable(comm_var)%vertical_layers

        ! go through the requested send pids
        ! and choose the requested
        DO np = 1, grid_comm_pattern%no_of_send_procs ! loop over PEs where to send the data
          IF (bfid /= grid_comm_pattern%send(np)%buffer_index) CYCLE

          DO var_no = 1, comm_variable(comm_var)%no_of_variables
            CALL fill_onepid_send_buffer_var(bfid, comm_var, np, var_no)
           END DO

        ENDDO

      END DO ! comm_var = 1, max_active_comm_variables

      buffer_start = send_procs_buffer(bfid)%start_index
      buffer_size  = send_procs_buffer(bfid)%buffer_size
      message_seq_id = 0
      DO WHILE (buffer_size > 0)

        message_size  = MIN(buffer_size, max_mpi_message_size)
        CALL p_isend(send_buffer(buffer_start:), &
          & send_procs_buffer(bfid)%pid, p_tag=halo_tag+message_seq_id, &
          &  p_count=message_size,              &
          & comm=my_work_communicator)

        buffer_size  = buffer_size  - message_size
        buffer_start = buffer_start + message_size
        message_seq_id = message_seq_id + 1

      ENDDO

    ENDDO !bfid = 1, active_send_buffers
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_fillandsend)

  END SUBROUTINE fill_and_nonblocksend_buffers
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Computes the sizes and the indexes for buffers to be sent
  SUBROUTINE compute_send_buffer_sizes()

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern

    INTEGER :: comm_var, vertical_layers, np, bfid, total_send_size, buffer_start
    INTEGER :: no_of_variables

    CHARACTER(*), PARAMETER :: method_name = "compute_send_buffer_sizes"

    ! clear buffer sizes
    DO bfid = 1, active_send_buffers
      send_procs_buffer(bfid)%buffer_size = 0
    ENDDO

    ! go through all variables and compute the requested send sizes
    DO comm_var = 1, max_active_comm_variables
      !
      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE

      grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
      vertical_layers = comm_variable(comm_var)%vertical_layers
      no_of_variables = comm_variable(comm_var)%no_of_variables

      ! go through the requested send pids
      DO np = 1, grid_comm_pattern%no_of_send_procs ! loop over PEs where to send the data

        total_send_size = &
          & ((grid_comm_pattern%send(np)%no_of_points * vertical_layers) &
          & + checksum_size) * no_of_variables

        bfid = grid_comm_pattern%send(np)%buffer_index
        send_procs_buffer(bfid)%buffer_size = send_procs_buffer(bfid)%buffer_size &
          &  + total_send_size

      ENDDO

    END DO

    ! now compute the buffers' start, current and end indexes
    buffer_start = 1
    DO bfid = 1, active_send_buffers
      send_procs_buffer(bfid)%start_index = buffer_start
      buffer_start = buffer_start + send_procs_buffer(bfid)%buffer_size
      send_procs_buffer(bfid)%end_index = buffer_start - 1
      ! set the current_index to start
      send_procs_buffer(bfid)%current_index = send_procs_buffer(bfid)%start_index
    ENDDO

    ! check if we went over the max_send_buffer_size
    IF (buffer_start >= max_send_buffer_size) &
      & CALL resize_send_buffer(buffer_start)
!       CALL finish(method_name, "buffer_start >= max_send_recv_buffer_size")

  END SUBROUTINE compute_send_buffer_sizes
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE fill_send_buffers_var(comm_var, var_no)
    INTEGER, INTENT(in) :: comm_var, var_no

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    REAL(wp), POINTER :: send_var_2d(:,:)
    REAL(wp), DIMENSION(:,:,:), POINTER :: send_var_3d

    INTEGER :: vertical_layers, np, bfid, var_send_size, current_buffer_index
    INTEGER :: i, k

!     CHARACTER(*), PARAMETER :: method_name = "fill_send_buffers_var"

    IF ( comm_variable(comm_var)%request /= communicate ) RETURN

    grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
    vertical_layers = comm_variable(comm_var)%vertical_layers
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      send_var_2d => comm_variable(comm_var)%send_values_2d(:,:)
    ELSEIF ( comm_variable(comm_var)%dim_4 > 0 ) THEN
      send_var_3d => comm_variable(comm_var)%send_values_4d(:,:,:,var_no)
    ELSE
      send_var_3d => comm_variable(comm_var)%send_values_3d
    ENDIF

    ! go through the requested send pids
    DO np = 1, grid_comm_pattern%no_of_send_procs ! loop over PEs where to send the data

      bfid = grid_comm_pattern%send(np)%buffer_index
      current_buffer_index = send_procs_buffer(bfid)%current_index
      var_send_size = grid_comm_pattern%send(np)%no_of_points * vertical_layers
      ! fill the header
      ! this is the id of the variable and the size we sent
      send_buffer(current_buffer_index) = REAL(comm_var,wp)
      current_buffer_index = current_buffer_index + 1
      send_buffer(current_buffer_index) = REAL(var_send_size,wp)
      current_buffer_index = current_buffer_index + 1

      ! fill the buffer
      IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
        ! fill 2d
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
          send_buffer(current_buffer_index) = send_var_2d &
            & ( grid_comm_pattern%send(np)%index_no(i),   &
                grid_comm_pattern%send(np)%block_no(i) )

            current_buffer_index = current_buffer_index + 1
        ENDDO

        IF (icon_comm_debug) THEN
          DO i = 1, grid_comm_pattern%send(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
              & send_procs_buffer(bfid)%pid, ":", &
              & i, send_var_2d( grid_comm_pattern%send(np)%index_no(i), &
                    grid_comm_pattern%send(np)%block_no(i) )
          ENDDO
        ENDIF

      ELSE

        ! fill 3d
!         DO k = 1, vertical_layers
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
            send_buffer(current_buffer_index : (current_buffer_index + vertical_layers - 1)) = &
               send_var_3d ( grid_comm_pattern%send(np)%index_no(i), 1:vertical_layers, &
                  grid_comm_pattern%send(np)%block_no(i) )

            current_buffer_index = current_buffer_index + vertical_layers
!             current_buffer_index = current_buffer_index + 1
        ENDDO
!         ENDDO

        IF (icon_comm_debug) THEN
          k=2
          DO i = 1, grid_comm_pattern%send(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
              & send_procs_buffer(bfid)%pid, ":", &
              & i,k, send_var_3d( grid_comm_pattern%send(np)%index_no(i), k, &
                    grid_comm_pattern%send(np)%block_no(i) )
          ENDDO
        ENDIF

      ENDIF

      send_procs_buffer(bfid)%current_index = current_buffer_index

    ENDDO

  END SUBROUTINE fill_send_buffers_var
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE fill_onepid_send_buffer_var(bfid, comm_var, send_proc_index, var_no)
    INTEGER, INTENT(in) :: bfid, comm_var, send_proc_index, var_no

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    REAL(wp), POINTER :: send_var_2d(:,:)
    REAL(wp), DIMENSION(:,:,:), POINTER :: send_var_3d

    INTEGER :: vertical_layers, np, var_send_size, current_buffer_index
    INTEGER :: i, k

!     CHARACTER(*), PARAMETER :: method_name = "fill_send_buffers_var"

    grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
    vertical_layers = comm_variable(comm_var)%vertical_layers
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      send_var_2d => comm_variable(comm_var)%send_values_2d(:,:)
    ELSEIF ( comm_variable(comm_var)%dim_4 > 0 ) THEN
      send_var_3d => comm_variable(comm_var)%send_values_4d(:,:,:,var_no)
    ELSE
      send_var_3d => comm_variable(comm_var)%send_values_3d
    ENDIF

    ! go through the requested send pids
    np = send_proc_index

    current_buffer_index = send_procs_buffer(bfid)%current_index
    var_send_size = grid_comm_pattern%send(np)%no_of_points * vertical_layers
    ! fill the header
    ! this is the id of the variabkle and the size we sent
    send_buffer(current_buffer_index) = REAL(comm_var,wp)
    send_buffer(current_buffer_index + 1) = REAL(var_send_size,wp)
    current_buffer_index = current_buffer_index + 2

    ! fill the buffer
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      ! fill 2d
      DO i = 1, grid_comm_pattern%send(np)%no_of_points
        send_buffer(current_buffer_index) = send_var_2d &
          & ( grid_comm_pattern%send(np)%index_no(i),   &
              grid_comm_pattern%send(np)%block_no(i) )

          current_buffer_index = current_buffer_index + 1
      ENDDO

      IF (icon_comm_debug) THEN
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
          write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
            & send_procs_buffer(bfid)%pid, ":", &
            & i, send_var_2d( grid_comm_pattern%send(np)%index_no(i), &
                  grid_comm_pattern%send(np)%block_no(i) )
        ENDDO
      ENDIF

    ELSE

      ! fill 3d
      DO i = 1, grid_comm_pattern%send(np)%no_of_points
          send_buffer(current_buffer_index : (current_buffer_index + vertical_layers - 1)) = &
            send_var_3d ( grid_comm_pattern%send(np)%index_no(i), 1:vertical_layers, &
                grid_comm_pattern%send(np)%block_no(i) )

          current_buffer_index = current_buffer_index + vertical_layers
      ENDDO

      IF (icon_comm_debug) THEN
        k=2
        DO i = 1, grid_comm_pattern%send(np)%no_of_points
          write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " sent to ", &
            & send_procs_buffer(bfid)%pid, ":", &
            & i,k, send_var_3d( grid_comm_pattern%send(np)%index_no(i), k, &
                  grid_comm_pattern%send(np)%block_no(i) )
        ENDDO
      ENDIF

    ENDIF ! (comm_variable(comm_var)%grid_dim == grid_2D )

    send_procs_buffer(bfid)%current_index = current_buffer_index


  END SUBROUTINE fill_onepid_send_buffer_var
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE nonblocksent_all_data(  )

    INTEGER :: bfid, buffer_start, buffer_size, message_size, message_seq_id
!     CHARACTER(*), PARAMETER :: method_name = "sent_all_buffers"

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_isend)

    DO bfid = 1, active_send_buffers

      buffer_start = send_procs_buffer(bfid)%start_index
      buffer_size  = send_procs_buffer(bfid)%buffer_size
      message_seq_id = 0
      DO WHILE ( buffer_size > 0 )

        message_size = MIN(buffer_size, max_mpi_message_size)
        CALL p_isend(send_buffer(buffer_start:), send_procs_buffer(bfid)%pid, &
          & p_tag=halo_tag+message_seq_id, p_count=message_size, &
          & comm=my_work_communicator)
        buffer_size  = buffer_size  - message_size
        buffer_start = buffer_start + message_size
        message_seq_id = message_seq_id + 1

      ENDDO

    ENDDO

    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_isend)

  END SUBROUTINE nonblocksent_all_data
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE blocksent_all_data(  )

    INTEGER :: bfid, buffer_start, buffer_size, message_size, message_seq_id
!     CHARACTER(*), PARAMETER :: method_name = "sent_all_buffers"

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_send)

!ICON_OMP_PARALLEL IF(icon_comm_openmp)
!ICON_OMP_DO PRIVATE(bfid, buffer_start, buffer_size, message_seq_id, message_size)
    DO bfid = 1, active_send_buffers

      buffer_start = send_procs_buffer(bfid)%start_index
      buffer_size  = send_procs_buffer(bfid)%buffer_size
      message_seq_id = 0
      DO WHILE ( buffer_size > 0 )

        message_size = MIN(buffer_size, max_mpi_message_size)
        CALL p_send(send_buffer(buffer_start:), send_procs_buffer(bfid)%pid, &
          & p_tag=halo_tag+message_seq_id, p_count=message_size, &
          & comm=my_work_communicator)
        buffer_size  = buffer_size  - message_size
        buffer_start = buffer_start + message_size
        message_seq_id = message_seq_id + 1

      ENDDO

    ENDDO
!ICON_OMP_END_DO_NOWAIT
!ICON_OMP_END_PARALLEL

    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_send)

  END SUBROUTINE blocksent_all_data
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  SUBROUTINE nonblockrecv_all_data(  )

    INTEGER :: bfid, buffer_start, buffer_size, message_size, message_seq_id

    CALL compute_recv_buffer_sizes()

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_ircv)

    ! start receive
    DO bfid = 1, active_recv_buffers

      buffer_start = recv_procs_buffer(bfid)%start_index
      buffer_size  = recv_procs_buffer(bfid)%buffer_size
      message_seq_id = 0
      DO WHILE (buffer_size > 0 )

        message_size =  MIN(buffer_size,max_mpi_message_size)
        CALL p_irecv(recv_buffer(buffer_start:), recv_procs_buffer(bfid)%pid, &
          & p_tag=halo_tag+message_seq_id, p_count=message_size, &
          & comm=my_work_communicator)

          buffer_size  = buffer_size  - message_size
          buffer_start = buffer_start + message_size
          message_seq_id = message_seq_id + 1
      ENDDO

    ENDDO
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_ircv)

  END SUBROUTINE nonblockrecv_all_data
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! set the current recv index to start
  SUBROUTINE set_recv_current_index()
    INTEGER :: bfid
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%current_index = recv_procs_buffer(bfid)%start_index
    ENDDO
  END SUBROUTINE set_recv_current_index
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  !>
  SUBROUTINE fill_vars_from_recv_buffers()

    INTEGER :: comm_var, var_no, bfid

    IF (activate_sync_timers) CALL timer_start(timer_icon_comm_fillrecv)

    ! set the current recv index to start
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%current_index = recv_procs_buffer(bfid)%start_index
    ENDDO

    ! go through all active variables and fill from buffer
    DO comm_var = 1, max_active_comm_variables
      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE

        DO var_no = 1, comm_variable(comm_var)%no_of_variables
          CALL fill_var_from_recv_buffers(comm_var, var_no)
        END DO

    ENDDO
    IF (activate_sync_timers) CALL timer_stop(timer_icon_comm_fillrecv)

  END SUBROUTINE fill_vars_from_recv_buffers
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  SUBROUTINE fill_var_from_recv_buffers(comm_var, var_no)
    INTEGER, INTENT(in) :: comm_var, var_no

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern
    REAL(wp), POINTER :: recv_var_2d(:,:)
    REAL(wp), DIMENSION(:,:,:), POINTER :: recv_var_3d

    INTEGER :: vertical_layers, np, bfid, var_recv_size, current_buffer_index
    INTEGER :: recv_var, recv_size
    INTEGER :: i, k

    CHARACTER(*), PARAMETER :: method_name = "fill_var_from_recv_buffers"

    IF ( comm_variable(comm_var)%request /= communicate ) RETURN

    grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
    vertical_layers = comm_variable(comm_var)%vertical_layers
    IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
      recv_var_2d => comm_variable(comm_var)%recv_values_2d(:,:)
    ELSEIF ( comm_variable(comm_var)%dim_4 > 0 ) THEN
      recv_var_3d => comm_variable(comm_var)%recv_values_4d(:,:,:,var_no)
    ELSE
      recv_var_3d => comm_variable(comm_var)%recv_values_3d
    ENDIF

    ! go through the requested receive pids
    DO np = 1, grid_comm_pattern%no_of_recv_procs ! loop over PEs from where to receive the data

      bfid = grid_comm_pattern%recv(np)%buffer_index
      current_buffer_index = recv_procs_buffer(bfid)%current_index
      var_recv_size = grid_comm_pattern%recv(np)%no_of_points * vertical_layers
      ! fill the header
      ! this is the id of the variable and the size we sent
      recv_var = INT(recv_buffer(current_buffer_index))
      current_buffer_index = current_buffer_index + 1
      recv_size = INT(recv_buffer(current_buffer_index))
      current_buffer_index = current_buffer_index + 1
      IF ( recv_var /= comm_var ) &
        CALL finish(method_name, 'recv_var /= comm_var')
      IF ( recv_size /= var_recv_size ) &
        CALL finish(method_name, 'recv_size /= var_recv_size')

      ! fill from the buffer
      IF (comm_variable(comm_var)%grid_dim == grid_2D ) THEN
         ! fill 2d
        DO i = 1, grid_comm_pattern%recv(np)%no_of_points
          recv_var_2d &
            & ( grid_comm_pattern%recv(np)%index_no(i),    &
                grid_comm_pattern%recv(np)%block_no(i) ) = &
          recv_buffer(current_buffer_index)

          current_buffer_index = current_buffer_index + 1
        ENDDO

        IF (icon_comm_debug) THEN
          DO i = 1, grid_comm_pattern%recv(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " recv from ", &
              & recv_procs_buffer(bfid)%pid, ":", &
              & i, recv_var_2d( grid_comm_pattern%recv(np)%index_no(i), &
                    grid_comm_pattern%recv(np)%block_no(i) )
          ENDDO
        ENDIF


      ELSE
        ! fill 3d
!         DO k = 1, vertical_layers
        DO i = 1, grid_comm_pattern%recv(np)%no_of_points
            recv_var_3d( grid_comm_pattern%recv(np)%index_no(i), 1:vertical_layers, &
                  grid_comm_pattern%recv(np)%block_no(i) ) = &
            recv_buffer(current_buffer_index : (current_buffer_index + vertical_layers - 1))

            current_buffer_index = current_buffer_index + vertical_layers
!             current_buffer_index = current_buffer_index + 1
!           ENDDO
        ENDDO

        IF (icon_comm_debug) THEN
          k=2
          DO i = 1, grid_comm_pattern%recv(np)%no_of_points
            write(log_file_id,*) TRIM(comm_variable(comm_var)%name), " recv from ", &
              & recv_procs_buffer(bfid)%pid, ":", &
              & i,k, recv_var_3d( grid_comm_pattern%recv(np)%index_no(i), k, &
                    grid_comm_pattern%recv(np)%block_no(i) )
          ENDDO
        ENDIF

      ENDIF

      recv_procs_buffer(bfid)%current_index = current_buffer_index

    ENDDO

  END SUBROUTINE fill_var_from_recv_buffers
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !>
  !! Computes the sizes and the indexes for buffers to be received
  SUBROUTINE compute_recv_buffer_sizes()

    TYPE(t_grid_comm_pattern), POINTER :: grid_comm_pattern

    INTEGER :: comm_var, vertical_layers, np, bfid, total_recv_size, buffer_start
    INTEGER :: no_of_variables

    CHARACTER(*), PARAMETER :: method_name = "compute_recv_buffer_sizes"

    ! clear buffer sizes
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%buffer_size = 0
    ENDDO

    ! go through all variables and compute the requested sizes
    DO comm_var = 1, max_active_comm_variables

      IF ( comm_variable(comm_var)%request /= communicate ) CYCLE

      grid_comm_pattern => comm_variable(comm_var)%grid_comm_pattern
      vertical_layers = comm_variable(comm_var)%vertical_layers
      no_of_variables = comm_variable(comm_var)%no_of_variables

      ! go through the requested recieve pids
      DO np = 1, grid_comm_pattern%no_of_recv_procs ! loop over PEs from where to receive the data

        total_recv_size = &
          & ((grid_comm_pattern%recv(np)%no_of_points * vertical_layers) + checksum_size) &
          & * no_of_variables
        bfid = grid_comm_pattern%recv(np)%buffer_index
        recv_procs_buffer(bfid)%buffer_size = recv_procs_buffer(bfid)%buffer_size &
          &  + total_recv_size

      ENDDO

    END DO

    ! now compute the buffers' start, current and end indexes
    buffer_start = 1
    DO bfid = 1, active_recv_buffers
      recv_procs_buffer(bfid)%start_index = buffer_start
      buffer_start = buffer_start + recv_procs_buffer(bfid)%buffer_size
      recv_procs_buffer(bfid)%end_index = buffer_start - 1
      ! set the current_index to start
      recv_procs_buffer(bfid)%current_index = recv_procs_buffer(bfid)%start_index
    ENDDO

    ! check if we went over the max_recv_buffer_size
    IF (buffer_start >= max_recv_buffer_size) &
      & CALL resize_recv_buffer(buffer_start)
!       CALL finish(method_name, "buffer_start >= max_send_recv_buffer_size")

  END SUBROUTINE compute_recv_buffer_sizes
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
!   SUBROUTINE icon_comm_show(var,  comm_pattern)
!     TYPE(t_comm_pattern), INTENT(IN) :: comm_pattern
!     REAL(wp), POINTER, INTENT(IN)   :: var(:,:,:)
!
!
!     INTEGER :: dim_2, np, current_buffer_index
!     INTEGER :: i, k, endidx
!
!    IF(this_is_mpi_sequential) RETURN
!
!     dim_2=SIZE(var,2)
!     dim_2=1
!
!     DO np = 1, comm_pattern%np_send
!
!        endidx = comm_pattern%send_startidx(np) + &
!          & comm_pattern%send_count(np) - 1
!
!         DO i = comm_pattern%send_startidx(np), endidx
!           DO k = 1, dim_2
!
!             write(file_id(my_mpi_work_id),*) " is to send to ", comm_pattern%pelist_send(np), ":", &
!               i,k, var( comm_pattern%send_src_idx(i), k, &
!                      comm_pattern%send_src_blk(i) )
!           ENDDO
!         ENDDO
!
!     ENDDO
!
!     DO np = 1, comm_pattern%np_recv
!
!         endidx = comm_pattern%recv_startidx(np) + &
!           & comm_pattern%recv_count(np) - 1
!
!         DO i = comm_pattern%recv_startidx(np), endidx
!           DO k = 1, dim_2
!
!             write(file_id(my_mpi_work_id),*) " got from ", comm_pattern%pelist_recv(np), ":", &
!               i,k, var( comm_pattern%recv_dst_idx(i), k, &
!                      comm_pattern%recv_dst_blk(i) )
!           ENDDO
!         ENDDO
!
!     ENDDO
!
!     DO i = 1, comm_pattern%n_send
!         DO k = 1, dim_2
!           write(file_id(my_mpi_work_id),*) " send out :", i, k, &
!              var(comm_pattern%send_src_idx(i),k,comm_pattern%send_src_blk(i))
!        ENDDO
!     ENDDO
!
!     DO i = 1, comm_pattern%n_pnts
!         DO k = 1, dim_2
!           write(file_id(my_mpi_work_id),*) " received :", i, k, &
!            var(comm_pattern%recv_dst_idx(i),k,comm_pattern%recv_dst_blk(i))
!        ENDDO
!    ENDDO
!
!
!   END SUBROUTINE icon_comm_show
  !-----------------------------------------------------------------------

  !> This is a utility routine for the parallel range-searching
  !  algorithm with GNATs, the problem is described below in the
  !  subroutine "mpi_reduce_mindistance_pts".
  SUBROUTINE mintype_minfct(invec, inoutvec, len)
    INTEGER :: len
    TYPE(t_mpi_mintype) :: invec(len), inoutvec(len)
    ! local variables
    INTEGER :: i

    DO i=1,len
      IF (invec(i)%rdist < inoutvec(i)%rdist) THEN
        inoutvec(i) = invec(i)
      ELSE IF (invec(i)%rdist == inoutvec(i)%rdist) THEN
        IF (invec(i)%glb_index < inoutvec(i)%glb_index) THEN
          inoutvec(i) = invec(i)
        END IF
      END IF
    END DO
  END SUBROUTINE mintype_minfct


  !> This is a utility routine for the parallel range-searching
  !  algorithm with GNATs, the problem is described below.
  !
  !  This routine has been placed in the module "mo_icon_comm_lib",
  !  since it requires low-level MPI library calls for a user-defined
  !  reduction operation.
  !
  !  Problem statement:
  !  For a given point in lon-lat space we might end up with more than
  !  one "minimum distance" triangle when searching in parallel. Thus,
  !  we must reduce the list of "minimum triangles" over all working
  !  PEs. There may occur the situation that we have two or more
  !  triangles with identical distance (e.g. at the poles). Then we
  !  choose between these triangles based on the global triangle
  !  index. This is necessary to ensure a well-defined behaviour for a
  !  varying number of MPI tasks.
  !
  !  Initial implementation: F. Prill, DWD (2012-09-03)
  !
  SUBROUTINE mpi_reduce_mindistance_pts(in, total_dim, comm)
    TYPE(t_mpi_mintype),  INTENT(INOUT) :: in(:)
    INTEGER,                 INTENT(IN) :: total_dim
    INTEGER,                 INTENT(IN) :: comm
    ! local variables:
    INTEGER  :: mpi_type(2), mpi_block(2), min_type, ierr, mintype_op


  END SUBROUTINE mpi_reduce_mindistance_pts

END MODULE mo_icon_comm_lib
