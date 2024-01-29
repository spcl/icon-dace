#ifdef __PGI
!pgi$g opt=1
#endif
! This module provides all routines for dividing patches
! (including interpolation state) and setting up communication.
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
! Much of the code in this module is duplicated due to different
! synchronization MPI RMA methods, which perform differently on
! different platforms,
!   if your system has slow passive target synchronization
!   (mpi_win_lock/mpi_win_unlock), set the HAVE_SLOW_PASSIVE_TARGET_ONESIDED
!   preprocessor macro, e.g. by adding -DHAVE_SLOW_PASSIVE_TARGET_ONESIDED to
!   the ICON FFLAGS
!   otherwise (fast passive target one-sided communication), no action is
!   required
! Defining this macro results in global arrays to be constructed with
! active target synchronization and extra dist_mult_array_rma_sync
! calls need to be executed before requested remote data can be
! accessed.
MODULE mo_setup_subdivision
  !
  !-------------------------------------------------------------------------
  !    ProTeX FORTRAN source: Style 2
  !    modified for ICON project, DWD/MPI-M 2006
  !-------------------------------------------------------------------------
  !
  USE mo_sync,               ONLY: global_max, global_min, global_sum
  USE mo_kind,               ONLY: wp
  USE mo_util_string,        ONLY: int2string
  USE mo_util_file,          ONLY: get_filename_noext
  USE mo_impl_constants,     ONLY: min_rlcell, max_rlcell,  &
    & min_rledge, max_rledge, min_rlvert, max_rlvert, max_phys_dom,  &
    & min_rlcell_int, min_rledge_int, min_rlvert_int, max_hw
  USE mo_math_constants,     ONLY: pi
  USE mo_exception,          ONLY: finish, message, warning

  USE mo_run_config,         ONLY: msg_level
  USE mo_io_units,           ONLY: filename_max
  USE mo_model_domain,       ONLY: t_patch, p_patch_local_parent, t_pre_patch, &
       c_num_edges, c_parent, c_phys_id, c_neighbor, c_edge, &
       c_vertex, c_center, c_refin_ctrl, e_parent, e_cell, &
       e_refin_ctrl, v_cell, v_num_edges, v_vertex, v_refin_ctrl
  USE mo_decomposition_tools,ONLY: t_grid_domain_decomp_info, &
    &                              get_local_index, get_valid_local_index, &
    &                              set_inner_glb_index, set_outer_glb_index, &
    &                              uniform_partition
  USE mo_mpi,                ONLY: proc_split, p_max, num_work_procs, my_process_is_work
#ifndef NOMPI
  USE mo_mpi,                ONLY: MPI_COMM_NULL, p_int, &
       mpi_in_place, mpi_success, mpi_sum, mpi_lor, p_bool
# ifndef NO_MPI_CHOICE_ARG
  USE mpi,                   ONLY: MPI_Allreduce
# endif
#endif
  USE mo_mpi,                ONLY: p_comm_work, p_int, &
    & p_pe_work, p_n_work, my_process_is_mpi_parallel, p_alltoall, p_alltoallv

  USE mo_parallel_config,       ONLY:  nproma, ldiv_phys_dom, set_nproma, &
    & division_method, division_file_name, n_ghost_rows, ignore_nproma_use_nblocks_c, nblocks_c, &
    & div_geometric, ext_div_from_file, write_div_to_file, use_div_from_file, proc0_shift, &
    & nproma_sub, nblocks_sub, ignore_nproma_sub_use_nblocks_sub

  USE mo_communication,      ONLY: blk_no, idx_no, idx_1d
  USE mo_grid_config,         ONLY: n_dom, n_dom_start, patch_weight, l_limited_area
  USE mo_alloc_patches,ONLY: allocate_basic_patch, allocate_remaining_patch, &
                             deallocate_pre_patch
  USE mo_decomposition_tools, ONLY: divide_cells_by_location, t_cell_info, &
    & sort_cell_info_by_cell_number, partidx_of_elem_uniform_deco
  USE mo_math_utilities,      ONLY: fxp_lat, fxp_lon
  USE mo_master_control,      ONLY: my_process_is_oceanic, my_process_is_waves
#ifndef NOMPI
  USE mo_divide_cells_by_location_mpi, ONLY: divide_cells_by_location_mpi, &
       init_divide_cells_by_location_mpi
#endif
  USE mo_grid_config,         ONLY: use_dummy_cell_closure
  USE mo_util_sort,           ONLY: quicksort
  USE mo_dist_dir,            ONLY: dist_dir_setup
  USE ppm_distributed_array,  ONLY: dist_mult_array, global_array_desc, &
    &                               dist_mult_array_new, &
    &                               dist_mult_array_local_ptr, &
    &                               dist_mult_array_expose, &
    &                               dist_mult_array_delete, &
    &                               dist_mult_array_get, &
    &                               ppm_int
#ifdef __BLOCK_GET__
  USE ppm_distributed_array,  ONLY: dist_mult_array_get_i4_blk, &
    &                               dist_mult_init_blk_comm, dist_mult_do_blk_comm, &
    &                               dist_mult_end_blk_comm, &
    &                               proc_lst, ind_lst, recv_array, &
    &                               dist_mult_array_get_dp_blk, &
    &                               dist_mult_init_blk_comm_dp, dist_mult_do_blk_comm_dp, &
    &                               dist_mult_end_blk_comm_dp, &
    &                               proc_lst_dp, ind_lst_dp, recv_array_dp
  USE mo_mpi,                 ONLY: mpi_address_kind
#endif
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
  USE ppm_distributed_array,  ONLY: sync_mode_active_target, &
    &                               dist_mult_array_rma_sync
#endif
  USE ppm_extents,            ONLY: extent, extent_start, extent_end, extent_size
  USE mo_read_netcdf_distributed, ONLY: t_distrib_read_data, distrib_nf_open, &
    &                                   distrib_nf_close, distrib_read, &
    &                                   delete_distrib_read, setup_distrib_read
  USE ppm_distributed_array,  ONLY: dist_mult_array, global_array_desc
  USE mo_util_uuid_types,     ONLY: uuid_string_length
  USE mo_util_uuid,           ONLY: uuid_unparse
  USE mo_read_netcdf_broadcast_2, ONLY: netcdf_open_input, netcdf_close, &
    &                                   netcdf_read_att_int
  USE mo_fortran_tools,       ONLY: t_ptr_2d_int
  USE mo_netcdf_errhandler, ONLY: nf
  USE mo_netcdf

  IMPLICIT NONE

  PRIVATE

  !modules interface-------------------------------------------
  !subroutines
  PUBLIC :: decompose_domain

  TYPE nb_flag_list_elem
    INTEGER, ALLOCATABLE :: idx(:)
    INTEGER, ALLOCATABLE :: owner(:)
  END TYPE nb_flag_list_elem

  !> module name string
  CHARACTER(*), PARAMETER :: modname = 'mo_setup_subdivision'

CONTAINS

  !------------------------------------------------------------------
  !>
  !!  Divide patches and interpolation states.
  !!
  !!
  !!  @note Despite its name, this subroutine is also called in
  !!        sequential runs where it simply copies (and initializes)
  !!        the patch data structure.
  !!
  SUBROUTINE decompose_domain( p_patch, p_patch_pre )
#if defined(_CRAYFTN) && _RELEASE_MAJOR == 12 && _RELEASE_MINOR == 0 && _RELEASE_PATCHLEVEL == 2
!dir$ optimize(-O0)
#endif

    TYPE(t_patch), INTENT(INOUT), TARGET :: p_patch(n_dom_start:)
    TYPE(t_pre_patch), INTENT(INOUT), TARGET :: p_patch_pre(n_dom_start:)

    CHARACTER(*), PARAMETER :: routine = TRIM("mo_subdivision:decompose_domain")

    ! Local variables:
    INTEGER :: jg, jgp, jc, jgc, n, &
      &        n_procs_decomp
    INTEGER :: nprocs(p_patch_pre(1)%n_childdom)
    REAL(wp) :: weight(p_patch_pre(1)%n_childdom)
    TYPE(t_pre_patch), POINTER :: wrk_p_parent_patch_pre
    ! LOGICAL :: l_compute_grid
    INTEGER :: my_process_type, order_type_of_halos
    INTEGER, POINTER :: local_ptr(:)
    TYPE(dist_mult_array) :: dist_cell_owner, dist_cell_owner_p
    INTEGER :: jg_normal

    CALL message(routine, 'start of domain decomposition')
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    CALL message(routine, 'Is this taking a long time?&
         & Consider removing macro definition HAVE_SLOW_PASSIVE&
         &_TARGET_ONESIDED from compile time settings.')
#else
    CALL message(routine, 'Is this taking a long time?&
         & Consider adding macro definition HAVE_SLOW_PASSIVE_TARGET_ONESIDED&
         & to compile time settings.')
#endif
    ! This subroutine has 2 different operating modes:
    !
    !
    ! - If the subroutine is called from a parallel run,
    !   then the domain is split according to the number of processors
    !
    ! - When called from a single processor run, the domain is
    !   not split but simply copied (and the p_patch_local_parent is
    !   initialized properly).
    !
    n_procs_decomp = p_n_work - proc0_shift


    ! -----------------------------------------------------------------------------
    ! Check for processor splitting

    ! Check if the processor set should be split for patches of the 1st generation.
    ! This is the case if patch_weight > 0 for at least one of the root's childs.
    ! For clarity, we require that patch_weight=0 for any other patch

    ! Get weights for 1st generation patches
    proc_split = .FALSE.
    DO jc = 1, p_patch_pre(1)%n_childdom
      jgc = p_patch_pre(1)%child_id(jc)
      weight(jc) = patch_weight(jgc)
      IF(weight(jc) > 0._wp) proc_split = .TRUE.
    ENDDO

    ! Check if weights for other patches are 0 (for clarity only)
    IF(patch_weight(1) /= 0._wp) &
      CALL finish(routine,'Weight for root patch must be 0')
    DO jg = 2, n_dom
      jgp = p_patch_pre(jg)%parent_id
      IF(jgp /= 1 .AND. patch_weight(jg) > 0._wp) &
        CALL finish(routine,'Weight for higher level patch must be 0')
    ENDDO

    ! Initialization of indicator if a patch is owned by a PE (always true except in case of processor splitting)
    p_patch(:)%domain_is_owned = .FALSE.
    ! -----------------------------------------------------------------------------

    IF(proc_split) THEN

      ! -----------------------------------------------------------------------------
      ! CASE 1: DECOMPOSE WITH PROCESSOR SPLITTING
      ! -----------------------------------------------------------------------------

      IF(p_pe_work==0) WRITE(0,*) 'Splitting processor grid for first level patches'
      IF(p_pe_work==0) WRITE(0,'(a,10f12.3)') 'Weights for first level patches:',weight(:)

      ! In this case, the working processor set must be at least as big
      ! as the number of childs of the root patch
      IF(p_patch_pre(1)%n_childdom > n_procs_decomp) &
        CALL finish(routine,'Too few procs for processor grid splitting')

      ! Get the number of procs per patch according to weight(:).
      ! Every patch gets at least 1 proc (of course!):
      nprocs(:) = 1

      ! The remaining procs are divided among patches similar to
      ! the d'Hondt method for elections

      DO n = p_patch_pre(1)%n_childdom+1, n_procs_decomp
        jg = MAXLOC(weight(:)/REAL(nprocs(:)+1,wp),1)
        nprocs(jg) = nprocs(jg)+1
      ENDDO

      IF(p_pe_work==0) THEN
        WRITE(0,*) 'Processor splitting:'
        DO jc = 1, p_patch_pre(1)%n_childdom
          jgc =  p_patch_pre(1)%child_id(jc)
          WRITE(0,'(a,i0,a,f10.3,a,i0,a,i0)')                                              &
            &   'Patch ',jgc,' weight ',weight(jc),' gets ',nprocs(jc),' of ',n_procs_decomp
          IF (nprocs(jc) <= 1) &
            CALL finish(routine,'Processor splitting requires at least 2 PEs per patch')
        ENDDO
      ENDIF

      ! Set proc0, n_proc for all patches ...

      ! ... for the root patch and patch 0 if it exists
      ! proc0_shift is 0 unless specified otherwise in the namelist

      p_patch(n_dom_start:1)%n_proc = n_procs_decomp
      p_patch(n_dom_start:1)%proc0  = proc0_shift

      p_patch(n_dom_start:1)%domain_is_owned = .TRUE.

      ! ... for 1st generation childs

      n = proc0_shift
      DO jc = 1, p_patch_pre(1)%n_childdom
        jgc = p_patch_pre(1)%child_id(jc)
        p_patch(jgc)%proc0  = n
        p_patch(jgc)%n_proc = nprocs(jc)
        ! Determine which PEs will have owned grid points in a given nested domain;
        ! a gridpoint-less PE0 gets a special treatment in order to avoid complications in writing output to stderr
        IF (p_pe_work >= p_patch(jgc)%proc0 .AND. p_pe_work < p_patch(jgc)%proc0+nprocs(jc) .OR. &
            jc == 1 .AND. p_pe_work < proc0_shift) THEN
          p_patch(jgc)%domain_is_owned = .TRUE.
        ENDIF
        n = n + nprocs(jc)
      ENDDO

      ! ... for deeper level descendants

      DO jg = 2, n_dom

        jgp = p_patch_pre(jg)%parent_id

        IF(jgp /= 1) THEN
          p_patch(jg)%n_proc = p_patch(jgp)%n_proc
          p_patch(jg)%proc0  = p_patch(jgp)%proc0

          p_patch(jg)%domain_is_owned = p_patch(jgp)%domain_is_owned
        ENDIF

      ENDDO

    ELSE

      ! -----------------------------------------------------------------------------
      ! CASE 2: STANDARD DECOMPOSITION, NO SPLITTING
      ! -----------------------------------------------------------------------------
      ! No splitting, proc0, n_proc are identical for all patches

      IF(p_pe_work==0) WRITE(0,*) 'No splitting of processor grid'
      p_patch(:)%n_proc = n_procs_decomp
      p_patch(:)%proc0  = proc0_shift

      p_patch(n_dom_start:n_dom)%domain_is_owned = .TRUE.
    ENDIF

#ifdef NOMPI
    p_patch(:)%comm = 0
#else
    p_patch(:)%comm = MPI_COMM_NULL
#endif
    p_patch(:)%rank = -1

    ! -----------------------------------------------------------------------------
    ! Divide patches
    ! -----------------------------------------------------------------------------

    DO jg_normal = n_dom_start, n_dom
      ! Reorder index to start with a "normal" patch in case there is a reduced radiation grid
      IF (n_dom_start == 0) THEN
        IF (jg_normal == 0) THEN
          jg = 1
        ELSE IF(jg_normal == 1) THEN
          jg = 0
        ELSE
          jg = jg_normal
        ENDIF
      ELSE
        jg=jg_normal
      ENDIF

      jgp = p_patch_pre(jg)%parent_id

      IF(jg == n_dom_start) THEN
        NULLIFY(wrk_p_parent_patch_pre) ! Must be NULL for patch
      ELSE
        wrk_p_parent_patch_pre => p_patch_pre(jgp)
      ENDIF

      ! Here comes the actual domain decomposition:
      ! Every cells gets assigned an owner.

      IF (my_process_is_mpi_parallel()) THEN
        ! Set division method, divide_for_radiation is only used for patch 0
        ! (...  but not in limited area mode because the load balancing for radiation grids
        !       does not make sense there.)
        CALL divide_patch_cells(p_patch_pre(jg), jg, p_patch(jg)%n_proc,               &
          &                     p_patch(jg)%proc0, dist_cell_owner, dist_cell_owner_p, &
          &                     wrk_p_parent_patch_pre,                                &
          &                     divide_for_radiation = ((jg == 0) .AND. .NOT. l_limited_area))
      ELSE
        ! trivial "decomposition"
        CALL create_dist_cell_owner(dist_cell_owner, &
          &                         p_patch_pre(jg)%n_patch_cells_g, &
          &                         p_patch_pre(jg)%dist_array_comm, &
          &                         p_patch_pre(jg)%dist_array_pes)
        CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_ptr)
        local_ptr = 0
        IF (jg > n_dom_start) &
          CALL divide_parent_cells(p_patch_pre(jg), p_patch_pre(jgp), &
            &                      dist_cell_owner, dist_cell_owner_p)
      END IF

      ! Please note: Previously, for jg==0 no ghost rows were set.
      ! Currently, we need ghost rows for jg==0 also for dividing the int state and grf state
      ! Have still to check if int state/grf state is needed at all for jg==0,
      ! if this is not the case, the ghost rows can be dropped again.


      ! order_type_of_halos = 0 order for parent (l_compute_grid = false)
      !                       1 order root grid  (l_compute_grid = true)
      !                       2 all halos go to the end, for ocean

      IF (my_process_is_oceanic()) THEN
        order_type_of_halos = 2
      ELSE
        order_type_of_halos = 1
      ENDIF

      CALL divide_patch(p_patch(jg), p_patch_pre(jg), dist_cell_owner, &
        &               n_ghost_rows, order_type_of_halos, p_pe_work)

      IF(jg > n_dom_start) THEN

        order_type_of_halos = 0

        CALL divide_patch(p_patch_local_parent(jg), p_patch_pre(jgp), &
          &               dist_cell_owner_p, 1, order_type_of_halos,  p_pe_work)

        CALL dist_mult_array_delete(dist_cell_owner_p)
      ENDIF

      CALL dist_mult_array_delete(dist_cell_owner)

    ENDDO

    ! The patches may be discarded now
    DO jg = n_dom_start, n_dom
      CALL deallocate_pre_patch( p_patch_pre(jg) )
    ENDDO

  END SUBROUTINE decompose_domain

  SUBROUTINE init_dist_cell_owner(dist_cell_owner, cell_owner)
    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner
    INTEGER, OPTIONAL, INTENT(IN) :: cell_owner(:)

    INTEGER, POINTER :: local_ptr(:)

    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_ptr)
    IF (PRESENT(cell_owner)) THEN
      local_ptr(:) = cell_owner(LBOUND(local_ptr,1):UBOUND(local_ptr,1))
    ELSE
      local_ptr(:) = 0
    END IF

    CALL dist_mult_array_expose(dist_cell_owner)

  END SUBROUTINE init_dist_cell_owner

  SUBROUTINE create_dist_cell_owner(dist_cell_owner, n_cells_g, &
    &                               dist_array_comm, dist_array_pes)

    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner
    INTEGER, INTENT(IN) :: n_cells_g, dist_array_comm
    TYPE(extent), INTENT(IN) :: dist_array_pes

    TYPE(global_array_desc) :: dist_cell_owner_desc(1)
    TYPE(extent) :: local_chunk(1,1)
    INTEGER :: dist_array_pes_start, dist_array_pes_size

    dist_cell_owner_desc(1)%a_rank = 1
    dist_cell_owner_desc(1)%rect(1)%first = 1
    dist_cell_owner_desc(1)%rect(1)%size = n_cells_g
    dist_cell_owner_desc(1)%element_dt = ppm_int

    dist_array_pes_start = extent_start(dist_array_pes) - 1
    dist_array_pes_size = extent_size(dist_array_pes)

    ! compute uniform partition
    local_chunk(1,1) = uniform_partition(dist_cell_owner_desc(1)%rect(1), &
        &                                dist_array_pes_size, &
        &                                p_pe_work + 1 - dist_array_pes_start)

    dist_cell_owner = dist_mult_array_new(dist_cell_owner_desc, local_chunk, &
      &                                   dist_array_comm, &
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      &          sync_mode=sync_mode_active_target &
#else
      &          cache_size=MIN(10, CEILING(SQRT(REAL(dist_array_pes_size)))) &
#endif
      )
  END SUBROUTINE create_dist_cell_owner

  !-----------------------------------------------------------------------------
  !! Divides the cells of a patch (in wrk_p_patch_pre) for parallelization.
  !!
  !! Outputs the subdivsion in cell_owner(:) which must already be allocated
  !! to the correct size (wrk_p_patch_pre%n_patch_cells_g)
  !!
  !! If  wrk_p_parent_patch_pre is associated, this indicates that the patch has a parent
  !! which has consquences for subdivision (cell with the same parent must not
  !! get to different PEs).

  SUBROUTINE divide_patch_cells(wrk_p_patch_pre, patch_no, n_proc, proc0, &
       &                        dist_cell_owner, dist_cell_owner_p, &
       &                        wrk_p_parent_patch_pre, divide_for_radiation)

    TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre
    INTEGER, INTENT(IN)  :: patch_no !> The patch number,
                                     ! used to identify patch specific decomposition
    INTEGER, INTENT(IN)  :: n_proc !> Number of processors for split
    INTEGER, INTENT(IN)  :: proc0  !> First processor of patch
    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner !> Cell division
    !> Cell division for parent (only computed when wrk_p_parent_patch_pre
    !! associated)
    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner_p
    TYPE(t_pre_patch), POINTER :: wrk_p_parent_patch_pre
    ! Private flag if patch should be divided for radiation calculation
    LOGICAL, INTENT(in) :: divide_for_radiation

    INTEGER, POINTER :: local_owner_ptr(:)
    INTEGER, ALLOCATABLE :: flag_c(:)
    CHARACTER(LEN=uuid_string_length) :: p_grid_uuid, grid_uuid
    CHARACTER(LEN=filename_max) :: use_division_file_name ! if ext_div_from_file
    LOGICAL :: div_from_file, div_file_exists
    INTEGER :: num_proc_in_div_file, ncid, ret, i

    div_from_file = .FALSE.

    IF (use_div_from_file .OR. &
        (division_method(patch_no) == ext_div_from_file)) THEN

      ! Please note: Unfortunately we cannot use p_io for doing I/O,
      ! since this might be the test PE which is never calling this routine
      ! (this is the case in the actual setup).

      IF (division_file_name(patch_no) == "") THEN
        IF (ASSOCIATED(wrk_p_parent_patch_pre)) THEN
          CALL uuid_unparse(wrk_p_parent_patch_pre%grid_uuid, p_grid_uuid)
        ELSE
          p_grid_uuid = "-";
        END IF
        CALL uuid_unparse(wrk_p_patch_pre%grid_uuid, grid_uuid)
        use_division_file_name = &
          & TRIM(get_filename_noext(wrk_p_patch_pre%grid_filename)) // "_" //&
          & TRIM(p_grid_uuid) // "_" // TRIM(grid_uuid) // '_cell_owner.nc'
      ELSE
        use_division_file_name = division_file_name(patch_no)
      ENDIF

      ! check whether file is available
      INQUIRE(FILE=use_division_file_name, EXIST=div_file_exists)

      ! check whether number of procs matches
      IF (div_file_exists) THEN

        ncid = netcdf_open_input(use_division_file_name)

        num_proc_in_div_file = netcdf_read_att_int(ncid, "cell_owner", "nprocs")

        ret = netcdf_close(ncid)
      ELSE
        num_proc_in_div_file = 0
      END IF

      IF (division_method(patch_no) == ext_div_from_file) THEN
        IF (.NOT. div_file_exists) &
          CALL finish('divide_patch','decomposition file not found')
        IF (num_proc_in_div_file /= p_n_work) &
          CALL finish('divide_patch','wrong number of procs in decomposition file')
      END IF

      div_from_file = num_proc_in_div_file == p_n_work

      IF (div_from_file) THEN
        CALL read_netcdf_decomposition(use_division_file_name, dist_cell_owner, &
          &                            wrk_p_patch_pre%n_patch_cells_g, &
          &                            wrk_p_patch_pre%dist_array_comm, &
          &                            wrk_p_patch_pre%dist_array_pes)

        ! Assign the cell owners of the current patch to the parent cells
        ! for the construction of the local parent, set "-1" elsewhere.
        IF (ASSOCIATED(wrk_p_parent_patch_pre)) &
          CALL divide_parent_cells(wrk_p_patch_pre, wrk_p_parent_patch_pre, &
            &                      dist_cell_owner, dist_cell_owner_p)
      END IF

    END IF

    IF (.NOT. div_from_file) THEN

      IF(division_method(patch_no) == div_geometric) THEN
        IF(ASSOCIATED(wrk_p_parent_patch_pre)) THEN
          ALLOCATE(flag_c(wrk_p_parent_patch_pre%cells%local_chunk(1,1)%first:&
                          wrk_p_parent_patch_pre%cells%local_chunk(1,1)%first+&
        &                 wrk_p_parent_patch_pre%cells%local_chunk(1,1)%size-1))
        ELSE
          ALLOCATE(flag_c(wrk_p_patch_pre%cells%local_chunk(1,1)%first: &
                          wrk_p_patch_pre%cells%local_chunk(1,1)%first + &
        &                 wrk_p_patch_pre%cells%local_chunk(1,1)%size - 1))
        END IF
      ELSE IF(division_method(patch_no) == ext_div_from_file) THEN
        CALL finish('divide_patch','invalid division method')
      END IF

      ! Built-in subdivison methods
      IF(ASSOCIATED(wrk_p_parent_patch_pre)) THEN

        ! Cells with the same parent must not go to different PEs.
        ! Thus we have to divide in reality the subset of the parent cells
        ! which cover the actual patch and then assign the ownership
        ! of the cells of the actual patch according to the parent cells

        ! Determine the subset of the parent patch covering the actual patch
        ! by flagging the according cells
        CALL compute_flag_c(flag_c, wrk_p_patch_pre, wrk_p_parent_patch_pre)

        ! Divide subset of patch
        ! Receives the PE  numbers for every cell
        CALL divide_subset_geometric(flag_c, n_proc, wrk_p_parent_patch_pre, &
             dist_cell_owner_p, ASSOCIATED(wrk_p_parent_patch_pre), &
             divide_for_radiation = divide_for_radiation)

        CALL dist_mult_array_local_ptr(dist_cell_owner_p, 1, local_owner_ptr)
        DO i = LBOUND(local_owner_ptr, 1), UBOUND(local_owner_ptr, 1)
          local_owner_ptr(i) = MERGE(local_owner_ptr(i) + proc0, -1, &
               local_owner_ptr(i) >= 0)
        END DO
        CALL dist_mult_array_expose(dist_cell_owner_p)

        ! Owners of the cells of the parent patch are now in dist_cell_owner_p.
        ! Set owners in current patch from this

        CALL compute_dist_owner_from_dist_owner_p( &
          dist_cell_owner, wrk_p_patch_pre, wrk_p_parent_patch_pre, &
          dist_cell_owner_p)

      ELSE

        ! Set subset flags where the "subset" is the whole patch
        flag_c = 1

        ! No parent patch, simply divide current patch

        ! Divide complete patch
        CALL divide_subset_geometric(flag_c, n_proc, wrk_p_patch_pre, &
             dist_cell_owner, ASSOCIATED(wrk_p_parent_patch_pre), &
             divide_for_radiation = divide_for_radiation)

        CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_owner_ptr)
        DO i = LBOUND(local_owner_ptr, 1), UBOUND(local_owner_ptr, 1)
          local_owner_ptr(i) = MERGE(local_owner_ptr(i) + proc0, -1, &
               local_owner_ptr(i) >= 0)
        END DO

      ENDIF

      CALL dist_mult_array_expose(dist_cell_owner)

      IF (write_div_to_file) THEN

        IF (division_file_name(patch_no) == "") THEN
          IF (ASSOCIATED(wrk_p_parent_patch_pre)) THEN
            CALL uuid_unparse(wrk_p_parent_patch_pre%grid_uuid, p_grid_uuid)
          ELSE
            p_grid_uuid = "-";
          END IF
          CALL uuid_unparse(wrk_p_patch_pre%grid_uuid, grid_uuid)
          use_division_file_name = &
            & TRIM(get_filename_noext(wrk_p_patch_pre%grid_filename)) // "_" // &
            & TRIM(p_grid_uuid) // "_" // TRIM(grid_uuid) // '_cell_owner.nc'
        ELSE
          use_division_file_name = division_file_name(patch_no)
        ENDIF

        CALL write_netcdf_decomposition(use_division_file_name, dist_cell_owner, &
          &                             wrk_p_patch_pre%n_patch_cells_g)
      END IF
    ENDIF ! division_method

  CONTAINS

    SUBROUTINE compute_flag_c(flag_c, wrk_p_patch_pre, wrk_p_parent_patch_pre)

      INTEGER, ALLOCATABLE :: flag_c(:)
      TYPE(t_pre_patch), INTENT(IN) :: wrk_p_patch_pre
      TYPE(t_pre_patch), INTENT(IN) :: wrk_p_parent_patch_pre

      INTEGER, POINTER :: local_parent_ptr(:), local_phys_id_ptr(:)
      INTEGER, ALLOCATABLE :: parent_phys_id_map(:,:), num_send(:), num_recv(:), &
        &                     send_displs(:), recv_displs(:), flag_c_prep(:,:)
      INTEGER :: parent_part_rank, i, dist_array_pes_size
      TYPE(extent) :: global_extent

      flag_c = 0

      CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_phys_id, &
        &                            local_phys_id_ptr)
      CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_parent, &
        &                            local_parent_ptr)

      ALLOCATE(parent_phys_id_map(2,LBOUND(local_parent_ptr,1):&
        &                           UBOUND(local_parent_ptr,1)))

      parent_phys_id_map(1, :) = local_parent_ptr
      parent_phys_id_map(2, :) = local_phys_id_ptr

      CALL quicksort(parent_phys_id_map(1, :), parent_phys_id_map(2, :))

      global_extent%first = 1
      global_extent%size = wrk_p_parent_patch_pre%n_patch_cells_g

      dist_array_pes_size = extent_size(wrk_p_patch_pre%dist_array_pes)

      ALLOCATE(num_send(dist_array_pes_size), &
        &      num_recv(dist_array_pes_size), &
        &      send_displs(dist_array_pes_size), &
        &      recv_displs(dist_array_pes_size))

      num_send = 0
      DO i = LBOUND(local_parent_ptr, 1), UBOUND(local_parent_ptr, 1)
        parent_part_rank = &
          partidx_of_elem_uniform_deco(global_extent, dist_array_pes_size, &
          &                            local_parent_ptr(i))
        num_send(parent_part_rank) = num_send(parent_part_rank) + 1
      END DO

      CALL p_alltoall(num_send, num_recv, &
        &             wrk_p_parent_patch_pre%dist_array_comm)

      send_displs(1) = 0
      recv_displs(1) = 0
      DO i = 2, dist_array_pes_size
        send_displs(i) = send_displs(i-1) + num_send(i-1)
        recv_displs(i) = recv_displs(i-1) + num_recv(i-1)
      END DO
      ALLOCATE(flag_c_prep(2,recv_displs(dist_array_pes_size) + &
        &                    num_recv(dist_array_pes_size)))
      CALL p_alltoallv(parent_phys_id_map, num_send, send_displs, &
        &              flag_c_prep, num_recv, recv_displs, &
        &              wrk_p_parent_patch_pre%dist_array_comm)

      DO i = 1, SIZE(flag_c_prep, 2)
        flag_c(flag_c_prep(1,i)) = MAX(1, flag_c_prep(2,i))
      END DO
    END SUBROUTINE compute_flag_c

    SUBROUTINE compute_dist_owner_from_dist_owner_p(dist_cell_owner, &
      wrk_p_patch_pre, wrk_p_parent_patch_pre, dist_cell_owner_p)

      TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner
      TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre
      TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_parent_patch_pre
      TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner_p

      INTEGER, POINTER :: local_parent_ptr(:), local_owner_p_ptr(:), &
        &                 local_owner_ptr(:)
      INTEGER :: i, parent_part_rank, dist_array_pes_size
      INTEGER, ALLOCATABLE :: parent_cell_idx(:), parent_cell_idx_permutation(:), &
        &                     num_send(:), num_recv(:), send_displs(:), &
        &                     recv_displs(:), inquired_parent_cell_idx(:)
      TYPE(extent) :: global_extent

      CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_parent, &
        &                            local_parent_ptr)
      ALLOCATE(parent_cell_idx(LBOUND(local_parent_ptr,1): &
        &                      UBOUND(local_parent_ptr,1)), &
        &      parent_cell_idx_permutation(LBOUND(local_parent_ptr,1): &
        &                                  UBOUND(local_parent_ptr,1)))
      parent_cell_idx = local_parent_ptr
      parent_cell_idx_permutation = &
        (/(i,i=LBOUND(parent_cell_idx,1),UBOUND(parent_cell_idx,1))/)

      CALL quicksort(parent_cell_idx, parent_cell_idx_permutation)

      dist_array_pes_size = extent_size(wrk_p_patch_pre%dist_array_pes)

      global_extent%first = 1
      global_extent%size = wrk_p_parent_patch_pre%n_patch_cells_g

      ALLOCATE(num_send(dist_array_pes_size), &
        &      num_recv(dist_array_pes_size), &
        &      send_displs(dist_array_pes_size), &
        &      recv_displs(dist_array_pes_size))
      num_send = 0
      DO i = LBOUND(parent_cell_idx,1), UBOUND(parent_cell_idx,1)
        parent_part_rank = &
          partidx_of_elem_uniform_deco(global_extent, &
          &                            dist_array_pes_size, &
          &                            parent_cell_idx(i))
        num_send(parent_part_rank) = num_send(parent_part_rank) + 1
      END DO

      CALL p_alltoall(num_send, num_recv, &
        &             wrk_p_parent_patch_pre%dist_array_comm)

      send_displs(1) = 0
      recv_displs(1) = 0
      DO i = 2, dist_array_pes_size
        send_displs(i) = send_displs(i-1) + num_send(i-1)
        recv_displs(i) = recv_displs(i-1) + num_recv(i-1)
      END DO
      ALLOCATE( &
        &  inquired_parent_cell_idx(recv_displs(dist_array_pes_size) + &
        &  num_recv(dist_array_pes_size)))
      CALL p_alltoallv(parent_cell_idx, num_send, send_displs, &
        &              inquired_parent_cell_idx, num_recv, recv_displs, &
        &              wrk_p_parent_patch_pre%dist_array_comm)

      CALL dist_mult_array_local_ptr(dist_cell_owner_p, 1, local_owner_p_ptr)

      inquired_parent_cell_idx = local_owner_p_ptr(inquired_parent_cell_idx)

      CALL p_alltoallv(inquired_parent_cell_idx, num_recv, recv_displs, &
        &              parent_cell_idx, num_send, send_displs, &
        &              wrk_p_parent_patch_pre%dist_array_comm)

      CALL create_dist_cell_owner(dist_cell_owner, &
        &                         wrk_p_patch_pre%n_patch_cells_g, &
        &                         wrk_p_patch_pre%dist_array_comm, &
        &                         wrk_p_patch_pre%dist_array_pes)
      CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_owner_ptr)

      local_owner_ptr(parent_cell_idx_permutation) = parent_cell_idx

    END SUBROUTINE compute_dist_owner_from_dist_owner_p

  END SUBROUTINE divide_patch_cells

  !-------------------------------------------------------------------------------------------------
  !>
  !! Sets the owner for the division of the cells of the parent patch
  !! with the same subdivision as the child

  SUBROUTINE divide_parent_cells(p_patch_pre, p_parent_patch_pre, &
    &                            dist_cell_owner, dist_cell_owner_p)

    !> Patch for which the parent should be divided
    TYPE(t_pre_patch), INTENT(INOUT) :: p_patch_pre
    !> Parent patch
    TYPE(t_pre_patch), INTENT(INOUT) :: p_parent_patch_pre
    !> Ownership of cells in p_patch_pre
    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner
    !> Output: Cell division for parent. Must be allocated to n_patch_cells of!
    !! the global parent
    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner_p

    INTEGER :: child_chunk_first, child_chunk_last, child_chunk_size, i
    INTEGER :: n_cells_parent_g, parent_part_rank, dist_array_pes_size
    INTEGER :: accum_recv, accum_send
    TYPE(global_array_desc) :: dist_cell_owner_p_desc(1)
    TYPE(extent) :: parent_chunk(1,1)
    INTEGER, POINTER :: parent_chunk_ptr(:), child_chunk_ptr(:), &
      &                 child_chunk_parent_ptr(:)
    INTEGER, ALLOCATABLE :: num_parent_owner_send(:), num_parent_owner_recv(:),&
      &                     parent_owner_send_displ(:), &
      &                     parent_owner_recv_displ(:), &
      &                     parent_cell_owner_map(:,:), parent_chunk_prep(:,:)

    ! get local chunk
    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, child_chunk_ptr)
    child_chunk_first = LBOUND(child_chunk_ptr, 1)
    child_chunk_size = SIZE(child_chunk_ptr)
    child_chunk_last = child_chunk_first + child_chunk_size - 1
    CALL dist_mult_array_local_ptr(p_patch_pre%cells%dist, c_parent, &
      &                            child_chunk_parent_ptr)

    n_cells_parent_g = p_parent_patch_pre%n_patch_cells_g

    dist_cell_owner_p_desc(1)%a_rank = 1
    dist_cell_owner_p_desc(1)%rect(1)%first = 1
    dist_cell_owner_p_desc(1)%rect(1)%size = n_cells_parent_g
    dist_cell_owner_p_desc(1)%element_dt = ppm_int

    ! get parent cells global index for all child cells of child chunk and their
    ! respective owners (parent cells can occur more than once)
    ALLOCATE(parent_cell_owner_map(2, child_chunk_first:child_chunk_last))
    parent_cell_owner_map(1,:) = child_chunk_ptr
    parent_cell_owner_map(2,:) = child_chunk_parent_ptr

    ! sort both arrays according to parent cell global index
    CALL quicksort(parent_cell_owner_map(2,:), parent_cell_owner_map(1,:))

    dist_array_pes_size = extent_size(p_patch_pre%dist_array_pes)

    ! compute number of parent cells per process
    ALLOCATE(num_parent_owner_send(0:dist_array_pes_size-1), &
      &      num_parent_owner_recv(0:dist_array_pes_size-1), &
      &      parent_owner_send_displ(0:dist_array_pes_size-1), &
      &      parent_owner_recv_displ(0:dist_array_pes_size-1))
    num_parent_owner_send = 0
    DO i = child_chunk_first, child_chunk_last
      ! compute for each parent cell global index the respective process
      parent_part_rank = &
        partidx_of_elem_uniform_deco(dist_cell_owner_p_desc(1)%rect(1), &
        &                            dist_array_pes_size, &
        &                            parent_cell_owner_map(2,i))-1
      num_parent_owner_send(parent_part_rank) = &
        num_parent_owner_send(parent_part_rank) + 1
    END DO

    ! alltoall number of parent cells per process
    CALL p_alltoall(num_parent_owner_send, num_parent_owner_recv, &
      &             p_parent_patch_pre%dist_array_comm)

    ! alltoallv parent cell global index and owners
    accum_recv = 0
    accum_send = 0
    DO i = 0, dist_array_pes_size-1
      parent_owner_recv_displ(i) = accum_recv
      parent_owner_send_displ(i) = accum_send
      accum_recv = accum_recv + num_parent_owner_recv(i)
      accum_send = accum_send + num_parent_owner_send(i)
    END DO
    ALLOCATE( &
      parent_chunk_prep( &
      & 2, parent_owner_recv_displ(dist_array_pes_size-1) + &
      &    num_parent_owner_recv(dist_array_pes_size-1)))
    CALL p_alltoallv(parent_cell_owner_map, num_parent_owner_send, &
      &              parent_owner_send_displ, &
      &              parent_chunk_prep, num_parent_owner_recv, &
      &              parent_owner_recv_displ, &
      &              p_parent_patch_pre%dist_array_comm)

    ! sort according to parent cell global index
    CALL quicksort(parent_chunk_prep(2,:), parent_chunk_prep(1,:))

    ! check that each parent cell global index occurs either 0 or 4 times
    IF (MOD(SIZE(parent_chunk_prep, 2), 4) /= 0) &
      CALL finish('divide_parent_cells','Incomplete parent cell encountered')
    DO i = 1, SIZE(parent_chunk_prep, 2), 4
      IF (ANY((parent_chunk_prep(1,i) /= parent_chunk_prep(1,i+1:i+3)) .OR. &
        &     (parent_chunk_prep(2,i) /= parent_chunk_prep(2,i+1:i+3)))) &
        CALL finish('divide_parent_cells','Incomplete parent cell encountered')
    END DO

    ! set up dist_array
    parent_chunk = p_parent_patch_pre%cells%local_chunk
    dist_cell_owner_p = dist_mult_array_new(dist_cell_owner_p_desc, &
      &                                     parent_chunk, &
      &                                     p_parent_patch_pre%dist_array_comm, &
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      &          sync_mode=sync_mode_active_target &
#else
      &          cache_size=MIN(10, CEILING(SQRT(REAL(dist_array_pes_size)))) &
#endif
      )
    CALL dist_mult_array_local_ptr(dist_cell_owner_p, 1, parent_chunk_ptr)
    parent_chunk_ptr = -1
    DO i = 1, SIZE(parent_chunk_prep, 2), 4
      parent_chunk_ptr(parent_chunk_prep(2,i)) = parent_chunk_prep(1,i)
    END DO

    CALL dist_mult_array_expose(dist_cell_owner_p)

  END SUBROUTINE divide_parent_cells

  !-------------------------------------------------------------------------------------------------
  !! Divides a patch (in wrk_p_patch_pre) for parallelization.
  !!
  !! Parameters:
  !! cell_owner          owner PE numbers of the global cells
  !! n_boundary_rows     number of boundary rows to be added
  !! order_type_of_halos 1= a "normal" grid for prognostic computations is processed,
  !!                     0= a local parent grid is processed
  !!                       in the first case, a finer distinction between different halo
  !!                       cell/edge/vertex levels is made to optimize communication
  !!                     2=move all halos to the end (for ocean)
  !!
  !! On exit, the entries of wrk_p_patch are set.

  SUBROUTINE divide_patch(wrk_p_patch, wrk_p_patch_pre, dist_cell_owner, &
    &                     n_boundary_rows, order_type_of_halos, my_proc)

    TYPE(t_patch), INTENT(INOUT) :: wrk_p_patch ! output patch, designated as INOUT because
                                                ! a few attributes are already set
    TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre ! out is required for
                                                        ! the distributed arrays

    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner
    INTEGER, INTENT(IN) :: n_boundary_rows
    INTEGER, INTENT(IN) :: order_type_of_halos
    INTEGER, INTENT(IN) :: my_proc

    ! local variables
    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::divide_patch'
    INTEGER :: i, j, jl, jb, je, jv, jg, jc

#ifdef __PGIC__
    TYPE(nb_flag_list_elem), ALLOCATABLE :: flag2_c_list(:), &
      flag2_v_list(:), flag2_e_list(:)
    INTEGER, ALLOCATABLE :: n2_ilev_c(:), n2_ilev_v(:), n2_ilev_e(:)
#else
    TYPE(nb_flag_list_elem) :: flag2_c_list(0:2*n_boundary_rows), &
         flag2_v_list(0:n_boundary_rows+1), flag2_e_list(0:2*n_boundary_rows+1)
    INTEGER :: n2_ilev_c(0:2*n_boundary_rows), n2_ilev_v(0:n_boundary_rows+1), &
         n2_ilev_e(0:2*n_boundary_rows+1)
#endif
    INTEGER, ALLOCATABLE :: owned_cells(:), owned_edges(:), owned_verts(:)
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    INTEGER, ALLOCATABLE :: glb_cell_neighbor(:,:), glb_cell_edge(:,:), &
      &                     glb_cell_vertex(:,:), cell_parent_idx(:), edge_parent_idx(:)
#else
    INTEGER :: glb_cell_neighbor, glb_cell_edge, &
      &        glb_cell_vertex, cell_parent_idx, edge_parent_idx
#endif
#ifdef __BLOCK_GET__
    INTEGER :: ierr
    INTEGER(MPI_ADDRESS_KIND) :: lb, ex
    INTEGER :: sz
#endif

    IF (msg_level >= 10)  CALL message(routine, 'dividing patch')

    CALL compute_owned_cells()

#ifdef __PGIC__
    ALLOCATE(flag2_c_list(0:2*n_boundary_rows), &
      &      flag2_v_list(0:n_boundary_rows+1), &
      &      flag2_e_list(0:2*n_boundary_rows+1), &
      &      n2_ilev_c(0:2*n_boundary_rows), &
      &      n2_ilev_v(0:n_boundary_rows+1), &
      &      n2_ilev_e(0:2*n_boundary_rows+1))
#endif
    CALL compute_flag_lists(flag2_c_list, flag2_v_list, flag2_e_list, &
      &                     n2_ilev_c, n2_ilev_v, n2_ilev_e, &
      &                     n_boundary_rows, &
    ! removed arguments because NEC has some problems with them...
    !  &                     owned_edges, owned_verts, &
      &                     order_type_of_halos)

    DO i = 1, 2*n_boundary_rows
      CALL sort_flag_list_by_owner(flag2_c_list(i))
    END DO
    DO i = 1, 2*n_boundary_rows+1
      CALL sort_flag_list_by_owner(flag2_e_list(i))
    END DO

    !-----------------------------------------------------------------------------------------------
    ! Get the number of cells/edges/verts and other data for patch allocation
    !-----------------------------------------------------------------------------------------------

    CALL prepare_patch(wrk_p_patch_pre, wrk_p_patch, &
         SUM(n2_ilev_c(:)), SUM(n2_ilev_e(:)), SUM(n2_ilev_v(:)))
    

    !-----------------------------------------------------------------------------------------------
    ! Set the global ownership for cells, edges and verts (needed for boundary exchange).
    ! Please note that opposed to cells, the global owner for edges/verts is
    ! not unique at the boundaries of the divided patch.
    ! To minimize load imbalance at runtime, we set the owner to the PE with the highest processor number
    ! which is participating at this edge/vert if both PE numbers are even or odd, otherwise
    ! the lowest processor number is chosen
    !-----------------------------------------------------------------------------------------------
    CALL dist_dir_setup(wrk_p_patch%cells%decomp_info%owner_dist_dir, &
      &                 flag2_c_list(0)%idx(1:n2_ilev_c(0)), &
      &                 wrk_p_patch%n_patch_cells_g, &
      &                 p_comm_work, p_pe_work, p_n_work)
    CALL dist_dir_setup(wrk_p_patch%verts%decomp_info%owner_dist_dir, &
      &                 owned_verts, wrk_p_patch%n_patch_verts_g, &
      &                 p_comm_work, p_pe_work, p_n_work)
    CALL dist_dir_setup(wrk_p_patch%edges%decomp_info%owner_dist_dir, &
      &                 owned_edges, wrk_p_patch%n_patch_edges_g, &
      &                 p_comm_work, p_pe_work, p_n_work)

    !-----------------------------------------------------------------------------------------------
    ! Get the indices of local cells/edges/verts within patch and vice versa.

    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_cells, &
         n_patch_cve_g = wrk_p_patch_pre%n_patch_cells_g, &
         patch_id = wrk_p_patch_pre%id, &
         nblks = wrk_p_patch%nblks_c, &
         npromz = wrk_p_patch%npromz_c, &
         max_cell_connectivity = wrk_p_patch%cells%max_connectivity, &
         min_rlcve = min_rlcell, &
         min_rlcve_int = min_rlcell_int, &
         max_rlcve = max_rlcell, &
         max_ilev = 2*n_boundary_rows, &
         max_hw_cve = 2*max_hw, &
         flag2_list = flag2_c_list, &
         n2_ilev = n2_ilev_c, &
         decomp_info = wrk_p_patch%cells%decomp_info, &
         start_index = wrk_p_patch%cells%start_index, &
         start_block = wrk_p_patch%cells%start_block, &
         end_index = wrk_p_patch%cells%end_index, &
         end_block = wrk_p_patch%cells%end_block, &
         start_g = wrk_p_patch_pre%cells%start, &
         end_g = wrk_p_patch_pre%cells%end, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .TRUE., &
         refin_ctrl = wrk_p_patch_pre%cells%dist, &
         refin_ctrl_aidx = c_refin_ctrl, &
         refinement_predicate = refine_cells)
    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_edges, &
         n_patch_cve_g = wrk_p_patch_pre%n_patch_edges_g, &
         patch_id = wrk_p_patch_pre%id, &
         nblks = wrk_p_patch%nblks_e, &
         npromz = wrk_p_patch%npromz_e, &
         max_cell_connectivity = wrk_p_patch%cells%max_connectivity, &
         min_rlcve = min_rledge, &
         min_rlcve_int = min_rledge_int, &
         max_rlcve = max_rledge, &
         max_ilev = 2*n_boundary_rows+1, &
         max_hw_cve = 2*max_hw + 1, &
         flag2_list = flag2_e_list, &
         n2_ilev = n2_ilev_e, &
         decomp_info = wrk_p_patch%edges%decomp_info, &
         start_index = wrk_p_patch%edges%start_index, &
         start_block = wrk_p_patch%edges%start_block, &
         end_index = wrk_p_patch%edges%end_index, &
         end_block = wrk_p_patch%edges%end_block, &
         start_g = wrk_p_patch_pre%edges%start, &
         end_g = wrk_p_patch_pre%edges%end, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .FALSE., &
         refin_ctrl = wrk_p_patch_pre%edges%dist, &
         refin_ctrl_aidx = e_refin_ctrl, &
         refinement_predicate = refine_edges)
    !---------------------------------------------------------------------------------------
    CALL build_patch_start_end(n_patch_cve = wrk_p_patch%n_patch_verts, &
         n_patch_cve_g = wrk_p_patch_pre%n_patch_verts_g, &
         patch_id = wrk_p_patch_pre%id, &
         nblks = wrk_p_patch%nblks_v, &
         npromz = wrk_p_patch%npromz_v, &
         max_cell_connectivity = wrk_p_patch%cells%max_connectivity, &
         min_rlcve = min_rlvert, &
         min_rlcve_int = min_rlvert_int, &
         max_rlcve = max_rlvert, &
         max_ilev = n_boundary_rows+1, &
         max_hw_cve = max_hw + 1, &
         flag2_list = flag2_v_list, &
         n2_ilev = n2_ilev_v, &
         decomp_info = wrk_p_patch%verts%decomp_info, &
         start_index = wrk_p_patch%verts%start_index, &
         start_block = wrk_p_patch%verts%start_block, &
         end_index = wrk_p_patch%verts%end_index, &
         end_block = wrk_p_patch%verts%end_block, &
         start_g = wrk_p_patch_pre%verts%start, &
         end_g = wrk_p_patch_pre%verts%end, &
         order_type_of_halos = order_type_of_halos, &
         l_cell_correction = .FALSE., &
         refin_ctrl = wrk_p_patch_pre%verts%dist, &
         refin_ctrl_aidx = v_refin_ctrl, &
         refinement_predicate = refine_verts)

    ! Sanity checks: are there still elements of the index lists filled with dummy values?
    ! Note: the use of write(0,*) instead of CALL message is intended here
    !       because the debug output is wanted for all PEs where an error occurs
    IF (ANY(wrk_p_patch%cells%start_index(:)<0) .OR. ANY(wrk_p_patch%cells%start_block(:)<0) .OR.&
        ANY(wrk_p_patch%cells%end_index(:)  <0) .OR. ANY(wrk_p_patch%cells%end_block(:)  <0)) THEN
      DO i = min_rlcell, max_rlcell
        WRITE(0,'(a,2i5,2i4,4i7)') 'cells',my_proc,wrk_p_patch%id,i,     &
          wrk_p_patch%cells%start_block(i),wrk_p_patch%cells%start_index(i), &
          wrk_p_patch%cells%end_block(i),  wrk_p_patch%cells%end_index(i)
      ENDDO
      CALL finish('divide_patch','Error in cell start/end indices')
    ENDIF

    IF (ANY(wrk_p_patch%edges%start_index(:)<0) .OR. ANY(wrk_p_patch%edges%start_block(:)<0) .OR.&
        ANY(wrk_p_patch%edges%end_index(:)  <0) .OR. ANY(wrk_p_patch%edges%end_block(:)  <0)) THEN
      DO i = min_rledge, max_rledge
        write(0,'(a,2i5,2i4,4i7)') 'edges',my_proc,wrk_p_patch%id,i,     &
          wrk_p_patch%edges%start_block(i),wrk_p_patch%edges%start_index(i), &
          wrk_p_patch%edges%end_block(i),  wrk_p_patch%edges%end_index(i)
      ENDDO
      CALL finish('divide_patch','Error in edge start/end indices')
    ENDIF

    IF (ANY(wrk_p_patch%verts%start_index(:)<0) .OR. ANY(wrk_p_patch%verts%start_block(:)<0) .OR.&
        ANY(wrk_p_patch%verts%end_index(:)  <0) .OR. ANY(wrk_p_patch%verts%end_block(:)  <0)) THEN
      DO i = min_rlvert, max_rlvert
        write(0,'(a,2i5,2i4,4i7)') 'verts',my_proc,wrk_p_patch%id,i,     &
          wrk_p_patch%verts%start_block(i),wrk_p_patch%verts%start_index(i), &
          wrk_p_patch%verts%end_block(i),  wrk_p_patch%verts%end_index(i)
      ENDDO
    CALL finish('divide_patch','Error in vertex start/end indices')
    ENDIF

    ! Fill 'old' two-dimensional index fields until removing them is completed
    DO j = 1, MAX(1,wrk_p_patch%n_childdom)
      wrk_p_patch%cells%start_blk(:,j) = wrk_p_patch%cells%start_block(:)
      wrk_p_patch%cells%start_idx(:,j) = wrk_p_patch%cells%start_index(:)
      wrk_p_patch%cells%end_blk(:,j)   = wrk_p_patch%cells%end_block(:)
      wrk_p_patch%cells%end_idx(:,j)   = wrk_p_patch%cells%end_index(:)

      wrk_p_patch%edges%start_blk(:,j) = wrk_p_patch%edges%start_block(:)
      wrk_p_patch%edges%start_idx(:,j) = wrk_p_patch%edges%start_index(:)
      wrk_p_patch%edges%end_blk(:,j)   = wrk_p_patch%edges%end_block(:)
      wrk_p_patch%edges%end_idx(:,j)   = wrk_p_patch%edges%end_index(:)

      wrk_p_patch%verts%start_blk(:,j) = wrk_p_patch%verts%start_block(:)
      wrk_p_patch%verts%start_idx(:,j) = wrk_p_patch%verts%start_index(:)
      wrk_p_patch%verts%end_blk(:,j)   = wrk_p_patch%verts%end_block(:)
      wrk_p_patch%verts%end_idx(:,j)   = wrk_p_patch%verts%end_index(:)
    ENDDO

    !-----------------------------------------------------------------------------------------------
    ! Set arrays of divided patch
    !-----------------------------------------------------------------------------------------------

    ! decomp_domain is -1 for invalid locations (at the end of the last strip)

    wrk_p_patch%cells%decomp_info%decomp_domain = -1
    wrk_p_patch%edges%decomp_info%decomp_domain = -1
    wrk_p_patch%verts%decomp_info%decomp_domain = -1

    !---------------------------------------------------------------------------------------

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    ALLOCATE(glb_cell_neighbor(wrk_p_patch%n_patch_cells, &
      &                        wrk_p_patch%cells%max_connectivity), &
      &      glb_cell_edge(wrk_p_patch%n_patch_cells, &
      &                    wrk_p_patch%cells%max_connectivity), &
      &      glb_cell_vertex(wrk_p_patch%n_patch_cells, &
      &                      wrk_p_patch%cells%max_connectivity))

#ifdef __BLOCK_GET__
    call dist_mult_init_blk_comm(3, wrk_p_patch%n_patch_cells, wrk_p_patch%cells%max_connectivity, wrk_p_patch_pre%cells%dist)
#endif

    DO j = 1, wrk_p_patch%n_patch_cells
      jg = wrk_p_patch%cells%decomp_info%glb_index(j)

      DO i=1,wrk_p_patch%cells%max_connectivity
#ifndef __BLOCK_GET__
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_neighbor, &
          &                      (/jg,i/), glb_cell_neighbor(j,i))
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
          &                      (/jg,i/), glb_cell_edge(j,i))
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
          &                      (/jg,i/), glb_cell_vertex(j,i))
#else
        CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_neighbor, &
          &                      (/jg,i/), glb_cell_neighbor(j,i), proc_lst(1,j,i), ind_lst(1,j,i))
        CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_edge, &
          &                      (/jg,i/), glb_cell_edge    (j,i), proc_lst(2,j,i), ind_lst(2,j,i))
        CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_vertex, &
          &                      (/jg,i/), glb_cell_vertex  (j,i), proc_lst(3,j,i), ind_lst(3,j,i))
#endif
      END DO
    END DO

#ifdef __BLOCK_GET__
    call dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
    CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)

#ifdef __BLOCK_GET__
    DO j = 1, wrk_p_patch%n_patch_cells
      DO i=1,wrk_p_patch%cells%max_connectivity
        IF (proc_lst(1,j,i) /= -1) THEN
          glb_cell_neighbor(j,i) = recv_array(proc_lst(1,j,i))%lst(ind_lst(1,j,i))
        END IF
        IF (proc_lst(2,j,i) /= -1) THEN
          glb_cell_edge    (j,i) = recv_array(proc_lst(2,j,i))%lst(ind_lst(2,j,i))
        END IF
        IF (proc_lst(3,j,i) /= -1) THEN
          glb_cell_vertex  (j,i) = recv_array(proc_lst(3,j,i))%lst(ind_lst(3,j,i))
        END IF
      END DO
    END DO

    call dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
#endif

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%cells%decomp_info%glb_index(j)

      DO i=1,wrk_p_patch%cells%max_connectivity
        CALL get_local_idx(wrk_p_patch%cells%decomp_info, &
          &                glb_cell_neighbor(j,i), jc)

        wrk_p_patch%cells%neighbor_idx(jl,jb,i) = idx_no(jc)
        wrk_p_patch%cells%neighbor_blk(jl,jb,i) = blk_no(jc)

        CALL get_local_idx(wrk_p_patch%edges%decomp_info, glb_cell_edge(j,i), je)
        wrk_p_patch%cells%edge_idx(jl,jb,i) = idx_no(je)
        wrk_p_patch%cells%edge_blk(jl,jb,i) = blk_no(je)

        CALL get_local_idx(wrk_p_patch%verts%decomp_info, &
          &                glb_cell_vertex(j, i), jv)

        wrk_p_patch%cells%vertex_idx(jl,jb,i) = idx_no(jv)
        wrk_p_patch%cells%vertex_blk(jl,jb,i) = blk_no(jv)

      ENDDO
    ENDDO

    DEALLOCATE(glb_cell_neighbor, glb_cell_edge, glb_cell_vertex)
#else

    DO j = 1, wrk_p_patch%n_patch_cells
      jg = wrk_p_patch%cells%decomp_info%glb_index(j)

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      DO i=1,wrk_p_patch%cells%max_connectivity
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_neighbor, &
          &                      (/jg,i/), glb_cell_neighbor)
        CALL get_local_idx(wrk_p_patch%cells%decomp_info, &
          &                glb_cell_neighbor, jc)

        wrk_p_patch%cells%neighbor_idx(jl,jb,i) = idx_no(jc)
        wrk_p_patch%cells%neighbor_blk(jl,jb,i) = blk_no(jc)

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
          &                      (/jg,i/), glb_cell_edge)
        CALL get_local_idx(wrk_p_patch%edges%decomp_info, &
             glb_cell_edge, je)

        wrk_p_patch%cells%edge_idx(jl,jb,i) = idx_no(je)
        wrk_p_patch%cells%edge_blk(jl,jb,i) = blk_no(je)

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
          &                      (/jg,i/), glb_cell_vertex)
        CALL get_local_idx(wrk_p_patch%verts%decomp_info, &
          &                glb_cell_vertex, jv)

        wrk_p_patch%cells%vertex_idx(jl,jb,i) = idx_no(jv)
        wrk_p_patch%cells%vertex_blk(jl,jb,i) = blk_no(jv)

      ENDDO
    ENDDO

#endif


#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED

    ALLOCATE(cell_parent_idx(wrk_p_patch%n_patch_cells))

#ifdef __BLOCK_GET__
    call dist_mult_init_blk_comm(3, wrk_p_patch%n_patch_cells, 1, wrk_p_patch_pre%cells%dist)
    call dist_mult_init_blk_comm_dp(2, wrk_p_patch%n_patch_cells, 1, wrk_p_patch_pre%cells%dist)
#endif

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%cells%decomp_info%glb_index(j)

      ! parent and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.

#ifndef __BLOCK_GET__
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_parent, (/jg/), &
        &                      cell_parent_idx(j))

      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, (/jg/), &
        &                      wrk_p_patch%cells%num_edges(jl,jb))
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, (/jg, 1/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lat)
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, (/jg, 2/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lon)
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_refin_ctrl, &
           (/jg/), wrk_p_patch%cells%refin_ctrl(jl,jb))
#else
      CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_parent, (/jg/), &
        &                      cell_parent_idx(j), proc_lst(1,j,1), ind_lst(1,j,1))
      CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_num_edges, (/jg/), &
        &                      wrk_p_patch%cells%num_edges(jl,jb), proc_lst(2,j,1), ind_lst(2,j,1))

      CALL dist_mult_array_get_dp_blk(wrk_p_patch_pre%cells%dist, c_center, (/jg, 1/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lat, proc_lst_dp(1,j,1), ind_lst_dp(1,j,1))
      CALL dist_mult_array_get_dp_blk(wrk_p_patch_pre%cells%dist, c_center, (/jg, 2/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lon, proc_lst_dp(2,j,1), ind_lst_dp(2,j,1))
      CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_refin_ctrl, &
           (/jg/), wrk_p_patch%cells%refin_ctrl(jl,jb), proc_lst(3,j,1), ind_lst(3,j,1))
#endif
    ENDDO

#ifdef __BLOCK_GET__
    CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
    CALL dist_mult_do_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif
    CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)

#ifdef __BLOCK_GET__
!$NEC ivdep
    DO j = 1, wrk_p_patch%n_patch_cells
      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (proc_lst(1,j,1) /= -1) THEN
        cell_parent_idx(j) = recv_array(proc_lst(1,j,1))%lst(ind_lst(1,j,1))
      END IF
      IF (proc_lst(2,j,1) /= -1) THEN
        wrk_p_patch%cells%num_edges(jl,jb) = recv_array(proc_lst(2,j,1))%lst(ind_lst(2,j,1))
      END IF
      IF (proc_lst_dp(1,j,1) /= -1) THEN
        wrk_p_patch%cells%center(jl,jb)%lat = recv_array_dp(proc_lst_dp(1,j,1))%lst(ind_lst_dp(1,j,1))
      END IF
      IF (proc_lst_dp(2,j,1) /= -1) THEN
        wrk_p_patch%cells%center(jl,jb)%lon = recv_array_dp(proc_lst_dp(2,j,1))%lst(ind_lst_dp(2,j,1))
      END IF
      IF (proc_lst(3,j,1) /= -1) THEN
        wrk_p_patch%cells%refin_ctrl(jl,jb) = recv_array(proc_lst(3,j,1))%lst(ind_lst(3,j,1))
      END IF
    END DO

    CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
    CALL dist_mult_end_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif

    DO j = 1, wrk_p_patch%n_patch_cells
      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch
      wrk_p_patch%cells%parent_glb_idx(jl,jb)  = idx_no(cell_parent_idx(j))
      wrk_p_patch%cells%parent_glb_blk(jl,jb)  = blk_no(cell_parent_idx(j))
    END DO

    DEALLOCATE(cell_parent_idx)

#else

    DO j = 1, wrk_p_patch%n_patch_cells

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%cells%decomp_info%glb_index(j)

      ! parent and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.

      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_parent, (/jg/), &
        &                      cell_parent_idx)
      wrk_p_patch%cells%parent_glb_idx(jl,jb)  = idx_no(cell_parent_idx)
      wrk_p_patch%cells%parent_glb_blk(jl,jb)  = blk_no(cell_parent_idx)

      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, (/jg/), &
        &                      wrk_p_patch%cells%num_edges(jl,jb))
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, (/jg, 1/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lat)
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, (/jg, 2/), &
        &                      wrk_p_patch%cells%center(jl,jb)%lon)
      CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_refin_ctrl, &
           (/jg/), wrk_p_patch%cells%refin_ctrl(jl,jb))
    ENDDO

#endif

    DO j = 0, 2 * n_boundary_rows
      DO i = 1, n2_ilev_c(j)
        jg = flag2_c_list(j)%idx(i)
        jc = get_local_index(wrk_p_patch%cells%decomp_info%glb2loc_index, jg)
        jb = blk_no(jc)
        jl = idx_no(jc)
        wrk_p_patch%cells%decomp_info%decomp_domain(jl,jb) = j
      END DO
    END DO

    !---------------------------------------------------------------------------------------

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED

    ALLOCATE(edge_parent_idx(wrk_p_patch%n_patch_edges))

#ifdef __BLOCK_GET__
    call dist_mult_init_blk_comm(2, wrk_p_patch%n_patch_edges, 1, wrk_p_patch_pre%edges%dist)
#endif
    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%edges%decomp_info%glb_index(j)

#ifndef __BLOCK_GET__
      CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_parent, &
        &                      (/jg/), edge_parent_idx(j))

      CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_refin_ctrl, &
           (/jg/), wrk_p_patch%edges%refin_ctrl(jl,jb))
#else
      CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_parent, &
        &                      (/jg/), edge_parent_idx(j), proc_lst(1,j,1), ind_lst(1,j,1))

      CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_refin_ctrl, &
           (/jg/), wrk_p_patch%edges%refin_ctrl(jl,jb), proc_lst(2,j,1), ind_lst(2,j,1))
#endif
    END DO

#ifdef __BLOCK_GET__
    CALL dist_mult_do_blk_comm(wrk_p_patch_pre%edges%dist)
#endif
    CALL dist_mult_array_rma_sync(wrk_p_patch_pre%edges%dist)

#ifdef __BLOCK_GET__
    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (proc_lst(1,j,1) /= -1) THEN
        edge_parent_idx(j) = recv_array(proc_lst(1,j,1))%lst(ind_lst(1,j,1))
      END IF
      IF (proc_lst(2,j,1) /= -1) THEN
        wrk_p_patch%edges%refin_ctrl(jl,jb) = recv_array(proc_lst(2,j,1))%lst(ind_lst(2,j,1))
      END IF

    ENDDO

    CALL dist_mult_end_blk_comm(wrk_p_patch_pre%edges%dist)
#endif

    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      ! parent and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.

      wrk_p_patch%edges%parent_glb_idx(jl,jb)    = idx_no(edge_parent_idx(j))
      wrk_p_patch%edges%parent_glb_blk(jl,jb)    = blk_no(edge_parent_idx(j))
    ENDDO

    DEALLOCATE(edge_parent_idx)

#else

    DO j = 1,wrk_p_patch%n_patch_edges

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%edges%decomp_info%glb_index(j)

      CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_parent, &
        &                      (/jg/), edge_parent_idx)
      wrk_p_patch%edges%parent_glb_idx(jl,jb)    = idx_no(edge_parent_idx)
      wrk_p_patch%edges%parent_glb_blk(jl,jb)    = blk_no(edge_parent_idx)

      ! parent and child_idx/child_blk still point to the global values.
      ! This will be changed in set_parent_child_relations.
      CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_refin_ctrl, &
           (/jg/), wrk_p_patch%edges%refin_ctrl(jl,jb))
    END DO

#endif

    DO j = 0, 2 * n_boundary_rows + 1
      DO i = 1, n2_ilev_e(j)
        jg = flag2_e_list(j)%idx(i)
        je = get_local_index(wrk_p_patch%edges%decomp_info%glb2loc_index, jg)
        jb = blk_no(je)
        jl = idx_no(je)
        wrk_p_patch%edges%decomp_info%decomp_domain(jl,jb) = j
      END DO
    END DO

    !---------------------------------------------------------------------------------------

#ifdef __BLOCK_GET__
    CALL dist_mult_init_blk_comm(1, wrk_p_patch%n_patch_verts, 1, wrk_p_patch_pre%verts%dist)
    CALL dist_mult_init_blk_comm_dp(2, wrk_p_patch%n_patch_verts, 1, wrk_p_patch_pre%verts%dist)
#endif
    DO j = 1,wrk_p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      jg = wrk_p_patch%verts%decomp_info%glb_index(j)

#ifndef __BLOCK_GET__
      CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
           (/jg, 1/), wrk_p_patch%verts%vertex(jl,jb)%lat)
      CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
           (/jg, 2/), wrk_p_patch%verts%vertex(jl,jb)%lon)
      CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_refin_ctrl, &
           (/jg/), wrk_p_patch%verts%refin_ctrl(jl,jb))
#else
      CALL dist_mult_array_get_dp_blk(wrk_p_patch_pre%verts%dist, v_vertex, &
           (/jg, 1/), wrk_p_patch%verts%vertex(jl,jb)%lat, proc_lst_dp(1,j,1), ind_lst_dp(1,j,1))
      CALL dist_mult_array_get_dp_blk(wrk_p_patch_pre%verts%dist, v_vertex, &
           (/jg, 2/), wrk_p_patch%verts%vertex(jl,jb)%lon, proc_lst_dp(2,j,1), ind_lst_dp(2,j,1))
      CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%verts%dist, v_refin_ctrl, &
           (/jg/), wrk_p_patch%verts%refin_ctrl(jl,jb), proc_lst(1,j,1), ind_lst(1,j,1))
#endif
    ENDDO

#ifdef __BLOCK_GET__
    CALL dist_mult_do_blk_comm(wrk_p_patch_pre%verts%dist)
    CALL dist_mult_do_blk_comm_dp(wrk_p_patch_pre%verts%dist)
#endif
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    CALL dist_mult_array_rma_sync(wrk_p_patch_pre%verts%dist)
#endif

#ifdef __BLOCK_GET__
    DO j = 1,wrk_p_patch%n_patch_verts

      jb = blk_no(j) ! Block index in distributed patch
      jl = idx_no(j) ! Line  index in distributed patch

      IF (proc_lst_dp(1,j,1) /= -1) THEN
        wrk_p_patch%verts%vertex(jl,jb)%lat = recv_array_dp(proc_lst_dp(1,j,1))%lst(ind_lst_dp(1,j,1))
      END IF
      IF (proc_lst_dp(2,j,1) /= -1) THEN
        wrk_p_patch%verts%vertex(jl,jb)%lon = recv_array_dp(proc_lst_dp(2,j,1))%lst(ind_lst_dp(2,j,1))
      END IF
      IF (proc_lst(1,j,1) /= -1) THEN
        wrk_p_patch%verts%refin_ctrl(jl,jb) = recv_array(proc_lst(1,j,1))%lst(ind_lst(1,j,1))
      END IF

    ENDDO

    CALL dist_mult_end_blk_comm(wrk_p_patch_pre%verts%dist)
    CALL dist_mult_end_blk_comm_dp(wrk_p_patch_pre%verts%dist)
#endif

    DO j = 0, n_boundary_rows + 1
      DO i = 1, n2_ilev_v(j)
        jg = flag2_v_list(j)%idx(i)
        jv = get_local_index(wrk_p_patch%verts%decomp_info%glb2loc_index, jg)
        jb = blk_no(jv)
        jl = idx_no(jv)
        wrk_p_patch%verts%decomp_info%decomp_domain(jl,jb) = j
      END DO
    END DO

    ! ------------------------------------------------------------
    ! > The following block masks negative refin_ctrl flags. This
    ! > happens just to test the correctness of the algorithm which
    ! > sets the parent-side of the refin_ctrl flags automatically.

    IF (msg_level >= 10) THEN
      CALL message(modname//': divide_patch', "erase negative refin_ctrl flags in patch "//&
        &TRIM(int2string(wrk_p_patch%id)))
    END IF

    DO j=1,wrk_p_patch%n_patch_cells
      jc = idx_no(j) ; jb = blk_no(j)
      IF (wrk_p_patch%cells%refin_ctrl(jc,jb) < 0)  wrk_p_patch%cells%refin_ctrl(jc,jb) = 0
    END DO
    DO j=1,wrk_p_patch%n_patch_edges
      jc = idx_no(j) ; jb = blk_no(j)
      IF (wrk_p_patch%edges%refin_ctrl(jc,jb) < 0)  wrk_p_patch%edges%refin_ctrl(jc,jb) = 0
    END DO
    DO j=1,wrk_p_patch%n_patch_verts
      jc = idx_no(j) ; jb = blk_no(j)
      IF (wrk_p_patch%verts%refin_ctrl(jc,jb) < 0)  wrk_p_patch%verts%refin_ctrl(jc,jb) = 0
    END DO

    ! ------------------------------------------------------------

    DEALLOCATE(owned_cells)

  CONTAINS

    SUBROUTINE compute_owned_cells()

      INTEGER, POINTER :: local_chunk(:)
      INTEGER :: local_chunk_first, local_chunk_last, local_chunk_size, i
      INTEGER, ALLOCATABLE :: n_cells_owned_per_proc(:)
      ! n_cells_from_proc(i) = number of cell ownerships transferred
      ! from process i
      INTEGER, ALLOCATABLE :: n_cells_from_proc(:)
      INTEGER, ALLOCATABLE :: glb_index_per_proc(:)
      INTEGER, ALLOCATABLE :: sorted_local_chunk(:)
      INTEGER, ALLOCATABLE :: send_displs(:)
      INTEGER, ALLOCATABLE :: recv_displs(:)
      INTEGER :: dist_array_pes_start, dist_array_pes_end, dist_array_pes_size

      dist_array_pes_start = extent_start(wrk_p_patch_pre%dist_array_pes) - 1
      dist_array_pes_end = extent_end(wrk_p_patch_pre%dist_array_pes) - 1
      dist_array_pes_size = extent_size(wrk_p_patch_pre%dist_array_pes)

      ! get local chunk
      CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_chunk)
      local_chunk_first = LBOUND(local_chunk, 1)
      local_chunk_size = SIZE(local_chunk)
      local_chunk_last = local_chunk_first + local_chunk_size - 1

      ! for parent patches most of the entries are -1
      ALLOCATE(n_cells_owned_per_proc(-1:p_n_work-1), &
        &      n_cells_from_proc(dist_array_pes_start:dist_array_pes_end), &
        &      glb_index_per_proc(local_chunk_first:local_chunk_last), &
        &      sorted_local_chunk(local_chunk_first:local_chunk_last))

      ! compute n_cells_owned_per_proc
      n_cells_owned_per_proc = 0
      DO i = local_chunk_first, local_chunk_last
        n_cells_owned_per_proc(local_chunk(i)) = &
          n_cells_owned_per_proc(local_chunk(i)) + 1

        ! initialise glb_index_per_proc
        glb_index_per_proc(i) = i
      END DO

      ! mpi_alltoall n_cells_owned_per_proc
      CALL p_alltoall(n_cells_owned_per_proc(dist_array_pes_start: &
        &                                    dist_array_pes_size-1), &
        &             n_cells_from_proc, &
        &             wrk_p_patch_pre%dist_array_comm)

      ! compute glb_index_per_proc sort by owner (see SUBROUTINE quicksort)
      sorted_local_chunk = local_chunk
      CALL quicksort(sorted_local_chunk, glb_index_per_proc)

      ! mpi_alltoallv glb_index_per_proc -> owned_cells
      ALLOCATE(send_displs(0:p_n_work-1), &
        &      recv_displs(dist_array_pes_start:dist_array_pes_end))
      send_displs(0) = n_cells_owned_per_proc(-1)
      recv_displs(dist_array_pes_start) = 0
      DO i = 1, p_n_work-1
        send_displs(i) = send_displs(i-1) + n_cells_owned_per_proc(i-1)
      END DO
      DO i = dist_array_pes_start + 1, dist_array_pes_end
        recv_displs(i) = recv_displs(i-1) + n_cells_from_proc(i-1)
      END DO
      ALLOCATE(owned_cells(SUM(n_cells_from_proc)))
      CALL p_alltoallv(glb_index_per_proc, &
        &              n_cells_owned_per_proc(dist_array_pes_start: &
        &                                     dist_array_pes_size-1), &
        &              send_displs(dist_array_pes_start: &
        &                          dist_array_pes_size-1), &
        &              owned_cells, n_cells_from_proc, &
        &              recv_displs, wrk_p_patch_pre%dist_array_comm)

      ! sort owned_cells
      CALL quicksort(owned_cells)

    END SUBROUTINE compute_owned_cells

    SUBROUTINE prepare_patch(wrk_p_patch_pre, wrk_p_patch, &
         n_patch_cells, n_patch_edges, n_patch_verts)
      !> output patch, designated as INOUT because
      !! a few attributes are already set
      TYPE(t_patch), INTENT(inout) :: wrk_p_patch
      TYPE(t_pre_patch), INTENT(in) :: wrk_p_patch_pre
      INTEGER, INTENT(in) :: n_patch_cells, n_patch_edges, n_patch_verts
      INTEGER :: new_nproma
      CHARACTER(LEN=64) :: message_text
      CHARACTER(LEN=64) :: routine = modname//'::prepare_patch'
      real(wp) :: nproma_avg
      INTEGER :: nproma_max, nproma_min

      wrk_p_patch%n_patch_cells = n_patch_cells
      IF (ignore_nproma_use_nblocks_c .AND. wrk_p_patch_pre%id == 1 ) THEN
        ! set nproma for work procs with gridpoints
        IF (my_process_is_work()) THEN
          IF(n_patch_cells > 0) THEN
            new_nproma = (n_patch_cells-1) / nblocks_c + 1
            CALL set_nproma(new_nproma)
          ENDIF          

          IF(num_work_procs > 1) THEN
            nproma_max = global_max(new_nproma)
            ! set nproma for work procs with no gridpoints
            IF(n_patch_cells < 1) THEN
              new_nproma = nproma_max
              CALL set_nproma(new_nproma)
            ENDIF
            nproma_min = global_min(new_nproma)
            nproma_avg = REAL(global_sum(new_nproma))/num_work_procs
            CALL message(routine, 'Set new nproma for each MPI work process.')
            write(message_text,'(a,2i7,f10.2)') 'max/min/avg: ', nproma_max, nproma_min, nproma_avg
            CALL message('# New nproma values: ', message_text)
          ELSE
            write(message_text,'(i7)') new_nproma
            CALL message('New nproma: ', message_text)
          ENDIF
          
        ENDIF

      ENDIF

      !! Also set the secondary nproma_sub and the number of subblocks nblocks_sub
      IF( ignore_nproma_sub_use_nblocks_sub ) THEN
        nproma_sub = (nproma-1)/nblocks_sub+1
      ELSE
        IF (nproma_sub > nproma) THEN
          CALL warning(routine, 'secondary nproma: nproma_sub cannot be bigger than nproma, adjusted')
          nproma_sub = nproma
        END IF
        nblocks_sub = (nproma-1)/nproma_sub+1
      END IF
      write(message_text,'(i7)') nproma_sub
      CALL message('Secondary nproma (nproma_sub): ', message_text)
      !$ACC UPDATE DEVICE(nproma_sub) ASYNC(1)

      wrk_p_patch%n_patch_edges = n_patch_edges
      wrk_p_patch%n_patch_verts = n_patch_verts

      ! save the number of cells/edges/verts of the patch
      wrk_p_patch%n_patch_cells_g = wrk_p_patch_pre%n_patch_cells_g
      wrk_p_patch%n_patch_edges_g = wrk_p_patch_pre%n_patch_edges_g
      wrk_p_patch%n_patch_verts_g = wrk_p_patch_pre%n_patch_verts_g
      !
      ! calculate and save values for the blocking, these are needed for patch allocation.
      ! NB: Avoid the case nblks=0 for empty patches, this might cause troubles
      ! if a empty patch is used somewhere (and npromz gets wrong in the formulas below).
      !
      ! ... for the cells
      wrk_p_patch%nblks_c       = blk_no(wrk_p_patch%n_patch_cells)
      wrk_p_patch%npromz_c      = wrk_p_patch%n_patch_cells - (wrk_p_patch%nblks_c - 1)*nproma
      wrk_p_patch%alloc_cell_blocks = wrk_p_patch%nblks_c
      IF (use_dummy_cell_closure) THEN
        IF (wrk_p_patch%npromz_c == nproma) &
             wrk_p_patch%alloc_cell_blocks = wrk_p_patch%nblks_c + 1
      ENDIF

      ! ... for the edges
      wrk_p_patch%nblks_e       = blk_no(wrk_p_patch%n_patch_edges)
      wrk_p_patch%npromz_e      = wrk_p_patch%n_patch_edges - &
        &                         (wrk_p_patch%nblks_e - 1) * nproma

      ! ... for the vertices
      wrk_p_patch%nblks_v       = blk_no(wrk_p_patch%n_patch_verts)
      wrk_p_patch%npromz_v      = wrk_p_patch%n_patch_verts - (wrk_p_patch%nblks_v - 1)*nproma

      ! Also needed for patch allocation
      wrk_p_patch%max_childdom  = wrk_p_patch_pre%max_childdom

      ! Set other scalar members of patch here too ..
      wrk_p_patch%grid_filename      = wrk_p_patch_pre%grid_filename
      wrk_p_patch%grid_filename_grfinfo = wrk_p_patch_pre%grid_filename_grfinfo
      wrk_p_patch%level              = wrk_p_patch_pre%level
      wrk_p_patch%nest_level         = wrk_p_patch_pre%nest_level
      wrk_p_patch%id                 = wrk_p_patch_pre%id
      wrk_p_patch%cells%max_connectivity = wrk_p_patch_pre%cells%max_connectivity
      wrk_p_patch%verts%max_connectivity = wrk_p_patch_pre%verts%max_connectivity
      wrk_p_patch%parent_id          = wrk_p_patch_pre%parent_id
      wrk_p_patch%parent_child_index = wrk_p_patch_pre%parent_child_index
      wrk_p_patch%child_id(:)        = wrk_p_patch_pre%child_id(:)
      wrk_p_patch%child_id_list(:)   = wrk_p_patch_pre%child_id_list(:)
      wrk_p_patch%n_childdom         = wrk_p_patch_pre%n_childdom
      wrk_p_patch%n_chd_total        = wrk_p_patch_pre%n_chd_total
      wrk_p_patch%nlev               = wrk_p_patch_pre%nlev
      wrk_p_patch%nlevp1             = wrk_p_patch_pre%nlevp1
      wrk_p_patch%nshift             = wrk_p_patch_pre%nshift
      wrk_p_patch%nshift_total       = wrk_p_patch_pre%nshift_total
      wrk_p_patch%nshift_child       = wrk_p_patch_pre%nshift_child
      wrk_p_patch%grid_uuid          = wrk_p_patch_pre%grid_uuid

      !-----------------------------------------------------------------------------------------------
      ! Allocate the required data arrays in patch
      !-----------------------------------------------------------------------------------------------

      CALL allocate_basic_patch(wrk_p_patch)
      CALL allocate_remaining_patch(wrk_p_patch,2) ! 2 = only those needed for parallelization control


    END SUBROUTINE prepare_patch

    ! this routine assumes that both lists are sorted in ascending order
    SUBROUTINE remove_entries_from_ref_list(list, n, ref_list)
      INTEGER, INTENT(inout) :: n
      INTEGER, INTENT(inout) :: list(n)
      INTEGER, INTENT(in) :: ref_list(:)

      INTEGER :: i, j, k, ref_list_size

      ref_list_size = SIZE(ref_list(:))

      IF (n > 0 .AND. ref_list_size > 0) THEN
        ! initialise j
        DO j = 1, ref_list_size
          IF (ref_list(j) >= list(1)) &
            EXIT
        END DO
        IF (j <= ref_list_size) THEN
          k = 0
          DO i = 1, n-1
            IF (list(i) /= ref_list(j)) THEN
              k = k + 1
              list(k) = list(i)
            END IF
            DO j = j, ref_list_size
              IF (ref_list(j) >= list(i+1)) &
                EXIT
            END DO
            IF (j > ref_list_size) EXIT
          END DO
          IF (j > ref_list_size) THEN
            DO i = i+1, n
              k = k + 1
              list(k) = list(i)
            END DO
          ELSE IF (list(i) /= ref_list(j)) THEN
            k = k + 1
            list(k) = list(i)
          END IF
          n = k
        END IF
      END IF
    END SUBROUTINE remove_entries_from_ref_list

    ! this routine assumes that the list is sorted in ascending order
    SUBROUTINE remove_duplicated_entries(list, n)
      INTEGER, INTENT(inout) :: n
      INTEGER, INTENT(inout) :: list(n)

      INTEGER :: i, j

      IF (n == 0) THEN
        RETURN
      END IF

      i = 1

      DO j = 2, n
        IF (list(j-1) /= list(j)) THEN
          i = i + 1
          list(i) = list(j)
        END IF
      END DO

      n = i

    END SUBROUTINE remove_duplicated_entries

    SUBROUTINE compute_flag_lists(flag2_c_list, flag2_v_list, flag2_e_list, &
      &                           n2_ilev_c, n2_ilev_v, n2_ilev_e, &
      &                           n_boundary_rows, &
    ! removed arguments because NEC has some problems with them...
    !  &                           owned_edges, owned_verts, &
      &                           order_type_of_halos)
      INTEGER, INTENT(IN) :: n_boundary_rows
      TYPE(nb_flag_list_elem), INTENT(OUT) :: &
        flag2_c_list(0:2*n_boundary_rows), flag2_v_list(0:n_boundary_rows+1), &
        flag2_e_list(0:2*n_boundary_rows+1)
      INTEGER, INTENT(OUT) :: n2_ilev_c(0:2*n_boundary_rows), n2_ilev_v(0:n_boundary_rows+1), &
        &  n2_ilev_e(0:2*n_boundary_rows+1)
      ! removed arguments because NEC has some problems with them...
      ! INTEGER, ALLOCATABLE, INTENT(OUT) :: owned_edges(:), owned_verts(:)
      INTEGER, INTENT(IN) :: order_type_of_halos

      INTEGER :: n, i, ic, jv, jv_, je, ilev, jc, k
      LOGICAL, ALLOCATABLE :: pack_mask(:)
      INTEGER, ALLOCATABLE :: temp_cells(:), temp_vertices(:), &
                              temp_vertices_owner(:), temp_edges(:), &
                              temp_edges_owner(:), inner_edges(:), &
                              edge_cells(:)
      INTEGER :: n_inner_edges, n_temp_edges, n_temp_cells, n_temp_vertices
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      INTEGER, ALLOCATABLE :: temp_num_edges(:), cell_edge(:,:), &
        &                     edge_cell(:,:,:), edge_cell_owner(:,:,:), &
        &                     cell_vertex(:,:), vert_cell(:,:)
#else
      INTEGER :: temp_num_edges, cell_edge, &
        & edge_cell, edge_cell_, edge_cell_owner, edge_cell_owner_, &
        &                     cell_vertex, vert_cell
#endif
#ifdef __BLOCK_GET__
      INTEGER :: tmp_n_temp_edges
#endif
      !---------------------------------------------------------------------
      ! flag_c_list(-1)%idx empty dummy list
      ! flag_c_list(0)%idx  all cells jg where cell_owner(jg) == my_proc
      ! flag_c_list(j)%idx j in 1..n_boundary_cells all cells bordering
      !                    on cells in flag_c_list(j - 1)
      ! flag_e_list(-1)%idx empty dummy list
      ! flag_e_list(j)%idx all edges where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the edge
      ! flag_v_list(-1)%idx empty dummy list
      ! flag_v_list(j)%idx all vertices where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the vertex
      !------------------------------------------------------------------------
      ! The purpose of the second set of flags is to control moving the
      ! halo points to the end of the index vector for nearly all cells
      ! even if order_type_of_halos = 1
      !------------------------------------------------------------------------
      ! flag2_e_list(0)%idx all edges which are owned by the local process
      ! flag2_e_list(1)%idx all edges adjacent to cells in flag_c_list(0)
      !                     but owned by other tasks
      ! flag2_e_list(2*i)%idx where i > 0, edges in idx are ONLY adjacent
      !                       to cells in flag_c_list(i)%idx
      ! flag2_e_list(2*i+1)%idx where i > 0, edges adjacent to one cell in
      !                         flag_c_list(i) and either
      !                            one cell in flag_c_list(i+1)
      !                         or edge is an outer edge,
      !                            i.e. adjacent to only one cell
      ! flag2_c_list(0)%idx == flag_c_list(0)%idx
      ! flag2_c_list(2*j-1)%idx where j > 0, contains all cells from
      !                         flag_c_list(j) which are
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_c_list(2*j)%idx   where j > 0, contains all cells from
      !                         flag_c_list(j) which are NOT
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_v_list(0)%idx contains all vertices from flag_v_list(0)%idx
      !                         owned by the current task
      ! flag2_v_list(1)%idx contains all vertices from flag_v_list(0)%idx
      !                     NOT owned by the current task
      ! flag2_v_list(j)%idx for j > 1 == flag_v_list(j-1)%idx
      !------------------------------------------------------------------------

      !--------------------------------------------------------------------------
      ! compute flag2 arrays
      !--------------------------------------------------------------------------

      ALLOCATE(pack_mask(0))

      ! collect cells of level 0
      n2_ilev_c(0) = SIZE(owned_cells(:))

      IF (n2_ilev_c(0) == wrk_p_patch_pre%n_patch_cells_g) THEN

        CALL compute_flag_lists_short(flag2_c_list, flag2_v_list, flag2_e_list, &
          &                           n2_ilev_c, n2_ilev_v, n2_ilev_e, &
          &                           n_boundary_rows &
        ! removed arguments because NEC has some problems with them...
        !  &                           ,owned_edges, owned_verts &
                                      )

        RETURN
      END IF

      ALLOCATE(flag2_c_list(0)%idx(n2_ilev_c(0)), &
        &      flag2_c_list(0)%owner(n2_ilev_c(0)))

      flag2_c_list(0)%idx(:) = owned_cells(:)
      flag2_c_list(0)%owner(:) = my_proc

      n_inner_edges = 0
      n_temp_edges = 0
      n_temp_vertices = 0
      ALLOCATE(inner_edges(n2_ilev_c(0) * &
                           wrk_p_patch_pre%cells%max_connectivity), &
               temp_edges(n2_ilev_c(0) * &
                          wrk_p_patch_pre%cells%max_connectivity), &
               edge_cells(n2_ilev_c(0) * &
                          wrk_p_patch_pre%cells%max_connectivity), &
               temp_vertices(n2_ilev_c(0) * &
                             wrk_p_patch_pre%cells%max_connectivity))

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED

      ALLOCATE(temp_num_edges(n2_ilev_c(0)), &
        &      cell_edge(n2_ilev_c(0), wrk_p_patch_pre%cells%max_connectivity), &
        &      edge_cell(n2_ilev_c(0), wrk_p_patch_pre%cells%max_connectivity, 2), &
        &      edge_cell_owner(n2_ilev_c(0), wrk_p_patch_pre%cells%max_connectivity, 2))

#ifdef __BLOCK_GET__
      CALL dist_mult_init_blk_comm(1, n2_ilev_c(0), 1, wrk_p_patch_pre%cells%dist)
#endif
      DO ic = 1, n2_ilev_c(0)

#ifndef __BLOCK_GET__
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
          &                      (/flag2_c_list(0)%idx(ic)/), &
          &                      temp_num_edges(ic))
#else
        CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_num_edges, &
          &                      (/flag2_c_list(0)%idx(ic)/), &
          &                      temp_num_edges(ic), proc_lst(1,ic,1), ind_lst(1,ic,1))
#endif
      END DO

#ifdef __BLOCK_GET__
      CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)

#ifdef __BLOCK_GET__
      DO ic = 1, n2_ilev_c(0)
        IF (proc_lst(1,ic,1) /= -1) THEN
          temp_num_edges(ic) = recv_array(proc_lst(1,ic,1))%lst(ind_lst(1,ic,1))
        END IF
      END DO

      CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)

      CALL dist_mult_init_blk_comm(1, n2_ilev_c(0), maxval(temp_num_edges(1:n2_ilev_c(0))), wrk_p_patch_pre%cells%dist)
#endif
      DO ic = 1, n2_ilev_c(0)

        DO i = 1, temp_num_edges(ic)

#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
            &                      (/flag2_c_list(0)%idx(ic),i/), &
            &                      cell_edge(ic, i))
#else
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_edge, &
            &                      (/flag2_c_list(0)%idx(ic),i/), &
            &                      cell_edge(ic, i), proc_lst(1,ic,i), ind_lst(1,ic,i))
#endif
        END DO
      END DO

#ifdef __BLOCK_GET__
      CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)

#ifdef __BLOCK_GET__
      DO ic = 1, n2_ilev_c(0)
        DO i = 1, temp_num_edges(ic)
          IF (proc_lst(1,ic,i) /= -1) THEN
            cell_edge(ic, i) = recv_array(proc_lst(1,ic,i))%lst(ind_lst(1,ic,i))
          END IF
        END DO
      END DO

      CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
      CALL dist_mult_init_blk_comm(2, n2_ilev_c(0), maxval(temp_num_edges(1:n2_ilev_c(0))), wrk_p_patch_pre%edges%dist)
#endif
      DO ic = 1, n2_ilev_c(0)

        DO i = 1, temp_num_edges(ic)

#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
               (/cell_edge(ic,i),1/), edge_cell(ic, i, 1))
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
               (/cell_edge(ic,i),2/), edge_cell(ic, i, 2))
#else
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_cell, &
               (/cell_edge(ic,i),1/), edge_cell(ic, i, 1), proc_lst(1,ic,i), ind_lst(1,ic,i))
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_cell, &
               (/cell_edge(ic,i),2/), edge_cell(ic, i, 2), proc_lst(2,ic,i), ind_lst(2,ic,i))
#endif
        END DO
      END DO

#ifdef __BLOCK_GET__
      CALL dist_mult_do_blk_comm(wrk_p_patch_pre%edges%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%edges%dist)

#ifdef __BLOCK_GET__
      DO ic = 1, n2_ilev_c(0)
        DO i = 1, temp_num_edges(ic)
          IF (proc_lst(1,ic,i) /= -1) THEN
            edge_cell(ic, i, 1) = recv_array(proc_lst(1,ic,i))%lst(ind_lst(1,ic,i))
          END IF
          IF (proc_lst(2,ic,i) /= -1) THEN
            edge_cell(ic, i, 2) = recv_array(proc_lst(2,ic,i))%lst(ind_lst(2,ic,i))
          END IF
        END DO
      END DO

      CALL dist_mult_end_blk_comm(wrk_p_patch_pre%edges%dist)
      CALL dist_mult_init_blk_comm(2, n2_ilev_c(0), maxval(temp_num_edges(1:n2_ilev_c(0))), dist_cell_owner)
#endif
      ! collect inner and outer edges and vertices adjacent to cells of level 0
      DO ic = 1, n2_ilev_c(0)

        DO i = 1, temp_num_edges(ic)

          IF (.NOT.(edge_cell(ic, i, 1) <= 0 .OR. &
            & edge_cell(ic, i, 2) <= 0 .OR. &
            & edge_cell(ic, i, 1) == edge_cell(ic, i, 2))) THEN

#ifndef __BLOCK_GET__
            CALL dist_mult_array_get(dist_cell_owner, 1, &
              &                      (/edge_cell(ic, i, 1)/), &
              &                      edge_cell_owner(ic, i, 1))
            CALL dist_mult_array_get(dist_cell_owner, 1, &
              &                      (/edge_cell(ic, i, 2)/), &
              &                      edge_cell_owner(ic, i, 2))
#else
            CALL dist_mult_array_get_i4_blk(dist_cell_owner, 1, &
              &                      (/edge_cell(ic, i, 1)/), &
              &                      edge_cell_owner(ic, i, 1), proc_lst(1,ic,i), ind_lst(1,ic,i))
            CALL dist_mult_array_get_i4_blk(dist_cell_owner, 1, &
              &                      (/edge_cell(ic, i, 2)/), &
              &                      edge_cell_owner(ic, i, 2), proc_lst(2,ic,i), ind_lst(2,ic,i))
#endif
          END IF
        END DO
      END DO
#ifdef __BLOCK_GET__
      call dist_mult_do_blk_comm(dist_cell_owner)
#endif
      CALL dist_mult_array_rma_sync(dist_cell_owner)
#ifdef __BLOCK_GET__
      DO ic = 1, n2_ilev_c(0)
        DO i = 1, temp_num_edges(ic)
          IF (proc_lst(1,ic,i) /= -1) THEN
            edge_cell_owner(ic, i, 1) = recv_array(proc_lst(1,ic,i))%lst(ind_lst(1,ic,i))
          END IF
          IF (proc_lst(2,ic,i) /= -1) THEN
            edge_cell_owner(ic, i, 2) = recv_array(proc_lst(2,ic,i))%lst(ind_lst(2,ic,i))
          END IF
        END DO
      END DO

      CALL dist_mult_end_blk_comm(dist_cell_owner)
      CALL dist_mult_init_blk_comm(1, n2_ilev_c(0), maxval(temp_num_edges(1:n2_ilev_c(0))), wrk_p_patch_pre%cells%dist)
#endif
      DO ic = 1, n2_ilev_c(0)

        DO i = 1, temp_num_edges(ic)

          IF (edge_cell(ic, i, 1) <= 0 .OR. &
            & edge_cell(ic, i, 2) <= 0 .OR. &
            & edge_cell(ic, i, 1) == edge_cell(ic, i, 2)) THEN
            n_inner_edges = n_inner_edges + 1
            inner_edges(n_inner_edges) = cell_edge(ic,i)
          ELSE
            IF (edge_cell_owner(ic,i,1) == edge_cell_owner(ic,i,2)) THEN
              n_inner_edges = n_inner_edges + 1
              inner_edges(n_inner_edges) = cell_edge(ic,i)
            ELSE
              n_temp_edges = n_temp_edges + 1
              temp_edges(n_temp_edges) = cell_edge(ic,i)
              edge_cells(n_temp_edges) = flag2_c_list(0)%idx(ic)
            END IF
          END IF

          n_temp_vertices = n_temp_vertices + 1
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
            &                      (/flag2_c_list(0)%idx(ic),i/), &
            &                      temp_vertices(n_temp_vertices))
#else
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_vertex, &
            &                      (/flag2_c_list(0)%idx(ic),i/), &
            &                      temp_vertices(n_temp_vertices), proc_lst(1,ic,i), ind_lst(1,ic,i))
#endif
        END DO
      END DO

#ifdef __BLOCK_GET__
      CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)

#ifdef __BLOCK_GET__
      n_temp_vertices = 0
      DO ic = 1, n2_ilev_c(0)
        DO i = 1, temp_num_edges(ic)
          n_temp_vertices = n_temp_vertices + 1
          IF (proc_lst(1,ic,i) /= -1) THEN
            temp_vertices(n_temp_vertices) = recv_array(proc_lst(1,ic,i))%lst(ind_lst(1,ic,i))
          END IF
        END DO
      END DO

      call dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
      DEALLOCATE(temp_num_edges, cell_edge, edge_cell, edge_cell_owner)

#else
      ! collect inner and outer edges and vertices adjacent to cells of level 0
      DO ic = 1, n2_ilev_c(0)

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
          &                      (/flag2_c_list(0)%idx(ic)/), temp_num_edges)
        DO i = 1, temp_num_edges

          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
            &                      (/flag2_c_list(0)%idx(ic),i/), cell_edge)

          CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
               (/cell_edge,1/), edge_cell)
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
               (/cell_edge,2/), edge_cell_)

          IF (      edge_cell  <= 0 &
               .OR. edge_cell_ <= 0 &
               .OR. edge_cell  == edge_cell_) THEN
            n_inner_edges = n_inner_edges + 1
            inner_edges(n_inner_edges) = cell_edge
          ELSE

            CALL dist_mult_array_get(dist_cell_owner, 1, (/edge_cell/), &
                 edge_cell_owner)
            CALL dist_mult_array_get(dist_cell_owner, 1, (/edge_cell_/), &
                 edge_cell_owner_)

            IF (edge_cell_owner == edge_cell_owner_) THEN
              n_inner_edges = n_inner_edges + 1
              inner_edges(n_inner_edges) = cell_edge
            ELSE
              n_temp_edges = n_temp_edges + 1
              temp_edges(n_temp_edges) = cell_edge
              edge_cells(n_temp_edges) = flag2_c_list(0)%idx(ic)
            END IF
          END IF

          n_temp_vertices = n_temp_vertices + 1
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
            &                      (/flag2_c_list(0)%idx(ic),i/), &
            &                      temp_vertices(n_temp_vertices))
        END DO
      END DO

#endif

      ! remove duplicated inner edges
      IF (n_inner_edges .GT. 1) THEN
        CALL quicksort(inner_edges(1:n_inner_edges))
        CALL remove_duplicated_entries(inner_edges(1:n_inner_edges), n_inner_edges)
      ENDIF
      ! remove duplicated vertices
      IF (n_temp_vertices .GT. 1) THEN
        CALL quicksort(temp_vertices(1:n_temp_vertices))
        CALL remove_duplicated_entries(temp_vertices(1:n_temp_vertices), &
                                     n_temp_vertices)
      ENDIF

      ! MoHa Note: vertices that are only adjacent to inner edges do not
      ! necessarily also have to be inner vertices...there are special cases...

      IF (SIZE(pack_mask(:)) < n_temp_vertices) THEN
        DEALLOCATE(pack_mask)
        ALLOCATE(pack_mask(n_temp_vertices))
      END IF

      ALLOCATE(temp_vertices_owner(n_temp_vertices))
      CALL compute_vertex_owner(temp_vertices(1:n_temp_vertices), &
        &                       temp_vertices_owner(1:n_temp_vertices))
      pack_mask(1:n_temp_vertices) = temp_vertices_owner(1:n_temp_vertices) &
        &                            == my_proc

      n = COUNT(pack_mask(1:n_temp_vertices))

      ALLOCATE(owned_verts(n))
      owned_verts(:) = PACK(temp_vertices(1:n_temp_vertices), &
                            pack_mask(1:n_temp_vertices))

      ! generate flag2_v_list(0) and flag2_v_list(1)
      IF (order_type_of_halos == 0 .OR. order_type_of_halos == 2 ) THEN
        ALLOCATE(flag2_v_list(0)%idx(n_temp_vertices), &
          &      flag2_v_list(0)%owner(n_temp_vertices), &
          &      flag2_v_list(1)%idx(0), &
          &      flag2_v_list(1)%owner(0))
        n2_ilev_v(0) = n_temp_vertices
        flag2_v_list(0)%idx(:) = temp_vertices(1:n_temp_vertices)
        flag2_v_list(0)%owner(:) = temp_vertices_owner(1:n_temp_vertices)
        n2_ilev_v(1) = 0
      ELSE
        ALLOCATE(flag2_v_list(0)%idx(n), flag2_v_list(0)%owner(n), &
                 flag2_v_list(1)%idx(n_temp_vertices - n), &
                 flag2_v_list(1)%owner(n_temp_vertices - n))
        n2_ilev_v(0) = n
        flag2_v_list(0)%idx(:) = owned_verts(:)
        flag2_v_list(0)%owner(:) = &
          PACK(temp_vertices_owner(1:n_temp_vertices), &
            &  pack_mask(1:n_temp_vertices))
        n2_ilev_v(1) = n_temp_vertices - n
        flag2_v_list(1)%idx(:) = PACK(temp_vertices(1:n_temp_vertices), &
                                      .NOT. pack_mask(1:n_temp_vertices))
        flag2_v_list(1)%owner(:) = &
          PACK(temp_vertices_owner(1:n_temp_vertices), &
            &  .NOT. pack_mask(1:n_temp_vertices))
      END IF

      CALL quicksort(temp_edges(1:n_temp_edges), edge_cells(1:n_temp_edges))
      IF (SIZE(pack_mask(:)) < n_temp_edges) THEN
        DEALLOCATE(pack_mask)
        ALLOCATE(pack_mask(n_temp_edges))
      END IF
      ALLOCATE(temp_edges_owner(n_temp_edges))
      CALL compute_edge_owner(temp_edges(1:n_temp_edges), &
        &                     temp_edges_owner(1:n_temp_edges))
      pack_mask(1:n_temp_edges) = &
        .NOT. temp_edges_owner(1:n_temp_edges) == my_proc
      n = COUNT(pack_mask(1:n_temp_edges))
      ALLOCATE(owned_edges(n_inner_edges + n_temp_edges - n))
      owned_edges(1:n_inner_edges) = inner_edges(1:n_inner_edges)
      owned_edges(n_inner_edges+1:) = PACK(temp_edges(1:n_temp_edges), &
                                             .NOT. pack_mask(1:n_temp_edges))
      CALL quicksort(owned_edges(:))

      ! generate flag2_e_list(0) and flag2_e_list(1)
  !    IF (n_boundary_rows > 0 .AND. order_type_of_halos /= 0) THEN
      IF (n_boundary_rows > 0 .AND. order_type_of_halos == 1) THEN

        ALLOCATE(flag2_e_list(0)%idx(n_inner_edges + n_temp_edges - n), &
          &      flag2_e_list(0)%owner(n_inner_edges + n_temp_edges - n), &
          &      flag2_e_list(1)%idx(n), flag2_e_list(1)%owner(n))
        n2_ilev_e(0) = n_inner_edges + n_temp_edges - n
        n2_ilev_e(1) = n
        flag2_e_list(0)%idx(:) = owned_edges(:)
        flag2_e_list(0)%owner(:) = my_proc
        flag2_e_list(1)%idx(:) = PACK(temp_edges(1:n_temp_edges), &
                                        pack_mask(1:n_temp_edges))
        flag2_e_list(1)%owner(:) = PACK(temp_edges_owner(1:n_temp_edges), &
                                        pack_mask(1:n_temp_edges))
      ELSE

        ALLOCATE(flag2_e_list(0)%idx(n_inner_edges + n_temp_edges), &
          &      flag2_e_list(0)%owner(n_inner_edges + n_temp_edges), &
          &      flag2_e_list(1)%idx(0), flag2_e_list(1)%owner(0))
        n2_ilev_e(0) = n_inner_edges + n_temp_edges
        n2_ilev_e(1) = 0
        flag2_e_list(0)%idx(1:n_inner_edges) = inner_edges(1:n_inner_edges)
        flag2_e_list(0)%owner(1:n_inner_edges) = my_proc
        flag2_e_list(0)%idx(n_inner_edges+1:) = temp_edges(1:n_temp_edges)
        flag2_e_list(0)%owner(n_inner_edges+1:) = temp_edges_owner(1:n_temp_edges)
        CALL quicksort(flag2_e_list(0)%idx(:), flag2_e_list(0)%owner(:))
      END IF

      DEALLOCATE(inner_edges)
      ALLOCATE(temp_cells(n_temp_edges))

      DO ilev = 1, n_boundary_rows

        ! collect cells of level 2*ilev-1
        n_temp_cells = 0
        IF (SIZE(temp_cells(:)) < n_temp_edges) THEN

          DEALLOCATE(temp_cells)
          ALLOCATE(temp_cells(n_temp_edges))
        END IF
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED

#ifdef __BLOCK_GET__
        CALL dist_mult_init_blk_comm(2, n_temp_edges, 1, wrk_p_patch_pre%edges%dist)
#endif
        ALLOCATE(edge_cell(n_temp_edges, 1, 2))
        DO i = 1, n_temp_edges
          je = temp_edges(i)
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
               (/je,1/), edge_cell(i,1,1))
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
               (/je,2/), edge_cell(i,1,2))
#else
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_cell, &
               (/je,1/), edge_cell(i,1,1), proc_lst(1,i,1), ind_lst(1,i,1))
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_cell, &
               (/je,2/), edge_cell(i,1,2), proc_lst(2,i,1), ind_lst(2,i,1))
#endif
        ENDDO

#ifdef __BLOCK_GET__
        CALL dist_mult_do_blk_comm(wrk_p_patch_pre%edges%dist)
#endif
        CALL dist_mult_array_rma_sync(wrk_p_patch_pre%edges%dist)

#ifdef __BLOCK_GET__
        DO i = 1, n_temp_edges
          IF (proc_lst(1,i,1) /= -1) THEN
            edge_cell(i,1,1) = recv_array(proc_lst(1,i,1))%lst(ind_lst(1,i,1))
          END IF
          IF (proc_lst(2,i,1) /= -1) THEN
            edge_cell(i,1,2) = recv_array(proc_lst(2,i,1))%lst(ind_lst(2,i,1))
          END IF
        END DO

        CALL dist_mult_end_blk_comm(wrk_p_patch_pre%edges%dist)
#endif
        DO i = 1, n_temp_edges

          jc = edge_cell(i,1,1)
          IF (jc == edge_cells(i)) THEN
            jc = edge_cell(i,1,2)
          END IF

          IF (jc > 0 .AND. jc /= edge_cells(i)) THEN
            n_temp_cells = n_temp_cells + 1
            temp_cells(n_temp_cells) = jc
          END IF
        END DO
        DEALLOCATE(edge_cell)

#else

        DO i = 1, n_temp_edges

          je = temp_edges(i)
          CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
               (/je,1/), edge_cell)

          IF (edge_cell == edge_cells(i)) THEN
            CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
                 (/je,2/),edge_cell)
          END IF

          IF (edge_cell > 0 .AND. edge_cell /= edge_cells(i)) THEN
            n_temp_cells = n_temp_cells + 1
            temp_cells(n_temp_cells) = edge_cell
          END IF
        END DO

#endif
        !remove duplicated entries
        IF (n_temp_cells .GT. 1) THEN
          CALL quicksort(temp_cells(1:n_temp_cells))
          CALL remove_duplicated_entries(temp_cells(1:n_temp_cells), n_temp_cells)
        ENDIF
        ! store cells of level 2*ilev-1
        n2_ilev_c(2*ilev-1) = n_temp_cells
        ALLOCATE(flag2_c_list(2*ilev-1)%idx(n_temp_cells), &
          &      flag2_c_list(2*ilev-1)%owner(n_temp_cells))
        flag2_c_list(2*ilev-1)%idx(:) = temp_cells(1:n_temp_cells)
#ifdef __BLOCK_GET__
        call dist_mult_init_blk_comm(1, n_temp_cells, 1, dist_cell_owner)
#endif
        DO i = 1, n_temp_cells
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(dist_cell_owner, 1, (/temp_cells(i)/), &
            &                      flag2_c_list(2*ilev-1)%owner(i))
#else
          CALL dist_mult_array_get_i4_blk(dist_cell_owner, 1, (/temp_cells(i)/), &
            &                      flag2_c_list(2*ilev-1)%owner(i), proc_lst(1,i,1), ind_lst(1,i,1))
#endif
        END DO
#ifdef __BLOCK_GET__
        CALL dist_mult_do_blk_comm(dist_cell_owner)
#endif
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
        CALL dist_mult_array_rma_sync(dist_cell_owner)
#endif

#ifdef __BLOCK_GET__
        DO i = 1, n_temp_cells
          IF (proc_lst(1,i,1) /= -1) THEN
            flag2_c_list(2*ilev-1)%owner(i) = recv_array(proc_lst(1,i,1))%lst(ind_lst(1,i,1))
          END IF
        END DO

        CALL dist_mult_end_blk_comm(dist_cell_owner)
#endif
        ! collect all cells adjacent to the outer vertices
        IF (SIZE(temp_cells(:)) < n_temp_vertices * &
          & wrk_p_patch_pre%verts%max_connectivity) THEN

          DEALLOCATE(temp_cells)
          ALLOCATE(temp_cells(n_temp_vertices * &
            &      wrk_p_patch_pre%verts%max_connectivity))
        END IF
        n_temp_cells = 0
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
        ALLOCATE(temp_num_edges(n_temp_vertices), &
          &      vert_cell(n_temp_vertices, &
          &                wrk_p_patch_pre%verts%max_connectivity))
#ifdef __BLOCK_GET__
        CALL dist_mult_init_blk_comm(1, n_temp_vertices, 1, wrk_p_patch_pre%verts%dist)
#endif
        DO jv = 1, n_temp_vertices
          jv_ = temp_vertices(jv)
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_num_edges, &
            &                      (/jv_/), temp_num_edges(jv))
#else
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%verts%dist, v_num_edges, &
            &                      (/jv_/), temp_num_edges(jv), proc_lst(1,jv,1), ind_lst(1,jv,1))
#endif
        END DO
#ifdef __BLOCK_GET__
        CALL dist_mult_do_blk_comm(wrk_p_patch_pre%verts%dist)
#endif
        CALL dist_mult_array_rma_sync(wrk_p_patch_pre%verts%dist)

#ifdef __BLOCK_GET__
        DO jv = 1, n_temp_vertices
          IF (proc_lst(1,jv,1) /= -1) THEN
            temp_num_edges(jv) = recv_array(proc_lst(1,jv,1))%lst(ind_lst(1,jv,1))
          END IF
        END DO

        CALL dist_mult_end_blk_comm(wrk_p_patch_pre%verts%dist)

        CALL dist_mult_init_blk_comm(1, n_temp_vertices, maxval(temp_num_edges(1:n_temp_vertices)), wrk_p_patch_pre%verts%dist)
#endif
        DO jv = 1, n_temp_vertices
          jv_ = temp_vertices(jv)
          DO i = 1, temp_num_edges(jv)
#ifndef __BLOCK_GET__
            CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_cell, &
                 (/jv_,i/), vert_cell(jv,i))
#else
            CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%verts%dist, v_cell, &
                 (/jv_,i/), vert_cell(jv,i), proc_lst(1,jv,i), ind_lst(1,jv,i))
#endif
          END DO
        END DO
#ifdef __BLOCK_GET__
        CALL dist_mult_do_blk_comm(wrk_p_patch_pre%verts%dist)
#endif
        CALL dist_mult_array_rma_sync(wrk_p_patch_pre%verts%dist)
#ifdef __BLOCK_GET__
        DO jv = 1, n_temp_vertices
          DO i = 1, temp_num_edges(jv)
            IF (proc_lst(1,jv,i) /= -1) THEN
              vert_cell(jv,i) = recv_array(proc_lst(1,jv,i))%lst(ind_lst(1,jv,i))
            END IF
          END DO
        END DO

        CALL dist_mult_end_blk_comm(wrk_p_patch_pre%verts%dist)
#endif
        DO jv = 1, n_temp_vertices
          jv_ = temp_vertices(jv)
          DO i = 1, temp_num_edges(jv)
            IF (vert_cell(jv,i) > 0) THEN
              n_temp_cells = n_temp_cells + 1
              temp_cells(n_temp_cells) = vert_cell(jv,i)
            END IF
          END DO
        END DO
        DEALLOCATE(temp_num_edges, vert_cell)

#else

        DO jv = 1, n_temp_vertices
          jv_ = temp_vertices(jv)
          CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_num_edges, &
            &                      (/jv_/), temp_num_edges)
          DO i = 1, temp_num_edges
            CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_cell, &
                 (/jv_,i/), vert_cell)
            IF (vert_cell > 0) THEN
              n_temp_cells = n_temp_cells + 1
              temp_cells(n_temp_cells) = vert_cell
            END IF
          END DO
        END DO

#endif

        IF (n_temp_cells .GT. 1) THEN
          CALL quicksort(temp_cells(1:n_temp_cells))
          CALL remove_duplicated_entries(temp_cells(1:n_temp_cells), n_temp_cells)
          ! remove cells that are on level 2*ilev-2 and level 2*ilev-1
          DO k = -2, -1
            CALL remove_entries_from_ref_list(temp_cells(1:n_temp_cells), &
              n_temp_cells, &
              flag2_c_list(2 * ilev + k)%idx(1:n2_ilev_c(2 * ilev + k)))
          END DO
        ENDIF
        ! store cells of level 2*ilev
        n2_ilev_c(2*ilev) = n_temp_cells
        ALLOCATE(flag2_c_list(2*ilev)%idx(n_temp_cells), &
          &      flag2_c_list(2*ilev)%owner(n_temp_cells))
        flag2_c_list(2*ilev)%idx(:) = temp_cells(1:n_temp_cells)
#ifdef __BLOCK_GET__
        CALL dist_mult_init_blk_comm(1, n_temp_cells, 1, dist_cell_owner)
#endif
        DO i = 1, n_temp_cells
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(dist_cell_owner, 1, (/temp_cells(i)/), &
            &                      flag2_c_list(2*ilev)%owner(i))
#else
          CALL dist_mult_array_get_i4_blk(dist_cell_owner, 1, (/temp_cells(i)/), &
            &                      flag2_c_list(2*ilev)%owner(i), proc_lst(1,i,1), ind_lst(1,i,1))
#endif
        END DO

#ifdef __BLOCK_GET__
        CALL dist_mult_do_blk_comm(dist_cell_owner)
#endif
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
        CALL dist_mult_array_rma_sync(dist_cell_owner)
#endif

#ifdef __BLOCK_GET__
        DO i = 1, n_temp_cells
          IF (proc_lst(1,i,1) /= -1) THEN
            flag2_c_list(2*ilev)%owner(i) = recv_array(proc_lst(1,i,1))%lst(ind_lst(1,i,1))
          END IF
        END DO

        CALL dist_mult_end_blk_comm(dist_cell_owner)
#endif
        ! get all edges of cells of level 2*ilev and level 2*ilev-1
        IF (SIZE(temp_edges(:)) < (n2_ilev_c(2*ilev) + &
          &                        n2_ilev_c(2*ilev - 1)) * &
          &                       wrk_p_patch_pre%cells%max_connectivity) THEN

          DEALLOCATE(temp_edges)
          ALLOCATE(temp_edges((n2_ilev_c(2*ilev) + &
            &                  n2_ilev_c(2*ilev - 1)) * &
            &                 wrk_p_patch_pre%cells%max_connectivity))
        END IF
        IF (SIZE(edge_cells(:)) < (n2_ilev_c(2*ilev) + &
          &                        n2_ilev_c(2*ilev - 1)) * &
          &                       wrk_p_patch_pre%cells%max_connectivity) THEN

          DEALLOCATE(edge_cells)
          ALLOCATE(edge_cells((n2_ilev_c(2*ilev) + &
            &                  n2_ilev_c(2*ilev - 1)) * &
            &                 wrk_p_patch_pre%cells%max_connectivity))
        END IF
        n_temp_edges = 0
        DO k = -1, 0
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
          ALLOCATE(temp_num_edges(n2_ilev_c(2*ilev+k)))
#ifdef __BLOCK_GET__
          CALL dist_mult_init_blk_comm(1, n2_ilev_c(2*ilev+k), 1, wrk_p_patch_pre%cells%dist)
#endif
          DO ic = 1, n2_ilev_c(2*ilev+k)
            jc = flag2_c_list(2*ilev+k)%idx(ic)
#ifndef __BLOCK_GET__
            CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges(ic))
#else
            CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges(ic), proc_lst(1,ic,1), ind_lst(1,ic,1))
#endif
          END DO
#ifdef __BLOCK_GET__
          call dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
          CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)
#ifdef __BLOCK_GET__
          DO ic = 1, n2_ilev_c(2*ilev+k)
            IF (proc_lst(1,ic,1) /= -1) THEN
              temp_num_edges(ic) = recv_array(proc_lst(1,ic,1))%lst(ind_lst(1,ic,1))
            END IF
          END DO

          CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
#endif

#ifdef __BLOCK_GET__
          CALL dist_mult_init_blk_comm(1,MAXVAL(temp_num_edges), &
                                       n2_ilev_c(2*ilev+k),wrk_p_patch_pre%cells%dist)
          tmp_n_temp_edges = n_temp_edges
#endif
          DO ic = 1, n2_ilev_c(2*ilev+k)
            jc = flag2_c_list(2*ilev+k)%idx(ic)
            DO i = 1, temp_num_edges(ic)
              n_temp_edges = n_temp_edges + 1
#ifndef __BLOCK_GET__
              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
                   (/jc,i/), temp_edges(n_temp_edges))
#else
              CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_edge, &
                   (/jc,i/), temp_edges(n_temp_edges), &
                   proc_lst(1,i,ic), ind_lst(1,i,ic))
#endif
              edge_cells(n_temp_edges) = jc
            END DO
          END  DO
#ifdef __BLOCK_GET__
          CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
          CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)

#ifdef __BLOCK_GET__
          n_temp_edges = tmp_n_temp_edges
          DO ic = 1, n2_ilev_c(2*ilev+k)
            jc = flag2_c_list(2*ilev+k)%idx(ic)
            DO i = 1, temp_num_edges(ic)
              n_temp_edges = n_temp_edges + 1
              IF (proc_lst(1,i,ic) /= -1) THEN
                temp_edges(n_temp_edges) = recv_array(proc_lst(1,i,ic))%lst(ind_lst(1,i,ic))
              END IF
            END DO
          END DO
          CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
          DEALLOCATE(temp_num_edges)
#else
          DO ic = 1, n2_ilev_c(2*ilev+k)
            jc = flag2_c_list(2*ilev+k)%idx(ic)
            CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges)
            DO i = 1, temp_num_edges
              n_temp_edges = n_temp_edges + 1
              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_edge, &
                   (/jc,i/), temp_edges(n_temp_edges))
              edge_cells(n_temp_edges) = jc
            END DO
          END DO
#endif
        END DO
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
        CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)
#endif
        IF (n_temp_edges > 1) THEN
          CALL quicksort(temp_edges(1:n_temp_edges), edge_cells(1:n_temp_edges))
          ! remove all edges of level 2*ilev-2 and 2*ilev-1
          DO k = -2, -1
            CALL remove_entries_from_ref_list(temp_edges(1:n_temp_edges), &
              n_temp_edges, flag2_e_list(2*ilev+k)%idx(1:n2_ilev_e(2*ilev+k)))
          END DO
          ! collect inner edges
          IF (SIZE(pack_mask(:)) < n_temp_edges) THEN
            DEALLOCATE(pack_mask)
            ALLOCATE(pack_mask(n_temp_edges))
          END IF
          ! mask inner edges (occur twice)
          pack_mask(1:n_temp_edges-1) = (temp_edges(1:n_temp_edges-1) == &
                                         temp_edges(2:n_temp_edges)) .AND. &
                                        (temp_edges(1:n_temp_edges-1) > 0)
          pack_mask(n_temp_edges) = .FALSE.
          n = COUNT(pack_mask(1:n_temp_edges))
          ALLOCATE(flag2_e_list(2*ilev)%idx(n), flag2_e_list(2*ilev)%owner(n))
          flag2_e_list(2*ilev)%idx(:) = PACK(temp_edges(1:n_temp_edges), &
                                               pack_mask(1:n_temp_edges))
          CALL compute_edge_owner(flag2_e_list(2*ilev)%idx(:), &
            &                     flag2_e_list(2*ilev)%owner(:))
          n2_ilev_e(2*ilev) = n
        ELSE
          ALLOCATE(flag2_e_list(2*ilev)%idx(0), flag2_e_list(2*ilev)%owner(0))
          n2_ilev_e(2*ilev) = 0
          pack_mask(1:n_temp_edges) = .FALSE.
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
          ! this call is needed for synchronization purposes
          CALL compute_edge_owner(flag2_e_list(2*ilev)%idx, &
            &                     flag2_e_list(2*ilev)%owner)
#endif
        END IF
        ! collect outer edges
        IF (n_temp_edges > 0) THEN
          ! mark inner edges
          pack_mask(2:n_temp_edges) = pack_mask(2:n_temp_edges) .OR. &
                                      pack_mask(1:n_temp_edges-1)
          ! mark invalid edges
          pack_mask(1:n_temp_edges) = pack_mask(1:n_temp_edges) .OR. &
                                      (temp_edges(1:n_temp_edges) <= 0)
          n = n_temp_edges - COUNT(pack_mask(1:n_temp_edges))
        ELSE
          n = 0
        END IF
        ALLOCATE(flag2_e_list(2*ilev+1)%idx(n), flag2_e_list(2*ilev+1)%owner(n))
        flag2_e_list(2*ilev+1)%idx(:) = PACK(temp_edges(1:n_temp_edges), &
                                               .NOT. pack_mask(1:n_temp_edges))
        CALL compute_edge_owner(flag2_e_list(2*ilev+1)%idx(:), &
          &                     flag2_e_list(2*ilev+1)%owner(:))
        n2_ilev_e(2*ilev+1) = n
        temp_edges(1:n) = flag2_e_list(2*ilev+1)%idx(:)
        edge_cells(1:n) = PACK(edge_cells(1:n_temp_edges), &
                               .NOT. pack_mask(1:n_temp_edges))
        n_temp_edges = n

        IF (SIZE(temp_vertices(:)) < SUM(n2_ilev_c(2*ilev-1:2*ilev)) * &
          &                              wrk_p_patch_pre%cells%max_connectivity) THEN

          DEALLOCATE(temp_vertices)
          ALLOCATE(temp_vertices(SUM(n2_ilev_c(2*ilev-1:2*ilev)) * &
            &                        wrk_p_patch_pre%cells%max_connectivity))
        END IF
        ! collect all vertices of level ilev + 1
        n_temp_vertices = 0
        DO k = -1, 0
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED

          ALLOCATE(temp_num_edges(n2_ilev_c(2*ilev+k)), &
            &      cell_vertex(n2_ilev_c(2*ilev+k), &
            &                  wrk_p_patch_pre%cells%max_connectivity))
#ifdef __BLOCK_GET__
          CALL dist_mult_init_blk_comm(1, n2_ilev_c(2*ilev+k), 1, wrk_p_patch_pre%cells%dist)
#endif
          DO ic = 1, n2_ilev_c(2*ilev+k)

            jc = flag2_c_list(2*ilev+k)%idx(ic)

#ifndef __BLOCK_GET__
            CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges(ic))
#else
            CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges(ic), proc_lst(1,ic,1), ind_lst(1,ic,1))
#endif
          END DO
#ifdef __BLOCK_GET__
          CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
          CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)
#ifdef __BLOCK_GET__
          DO ic = 1, n2_ilev_c(2*ilev+k)
            IF (proc_lst(1,ic,1) /= -1) THEN
              temp_num_edges(ic) = recv_array(proc_lst(1,ic,1))%lst(ind_lst(1,ic,1))
            END IF
          END DO

          CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
          CALL dist_mult_init_blk_comm(1, n2_ilev_c(2*ilev+k), MAXVAL(temp_num_edges(1:n2_ilev_c(2*ilev+k))), &
                                       wrk_p_patch_pre%cells%dist)
#endif
          DO ic = 1, n2_ilev_c(2*ilev+k)

            jc = flag2_c_list(2*ilev+k)%idx(ic)

            DO i = 1, temp_num_edges(ic)
#ifndef __BLOCK_GET__

              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/jc,i/), cell_vertex(ic, i))
#else
              CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/jc,i/), cell_vertex(ic, i), proc_lst(1,ic,i), ind_lst(1,ic,i))
#endif
            END DO
          END DO
#ifdef __BLOCK_GET__
          CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
          CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)
#ifdef __BLOCK_GET__
          DO ic = 1, n2_ilev_c(2*ilev+k)
            DO i = 1, temp_num_edges(ic)
              IF (proc_lst(1,ic,i) /= -1) THEN
                cell_vertex(ic, i) = recv_array(proc_lst(1,ic,i))%lst(ind_lst(1,ic,i))
              END IF
            END DO
          END DO

          CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
          DO ic = 1, n2_ilev_c(2*ilev+k)
            DO i = 1, temp_num_edges(ic)
              IF (cell_vertex(ic,i) > 0) THEN
                n_temp_vertices = n_temp_vertices + 1
                temp_vertices(n_temp_vertices) = cell_vertex(ic, i)
              END IF
            END DO
          END DO
          DEALLOCATE(temp_num_edges, cell_vertex)

#else

          DO ic = 1, n2_ilev_c(2*ilev+k)

            jc = flag2_c_list(2*ilev+k)%idx(ic)

            CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_num_edges, &
              &                      (/jc/), temp_num_edges)

            DO i = 1, temp_num_edges

              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/jc,i/),cell_vertex)

              IF (cell_vertex > 0) THEN
                n_temp_vertices = n_temp_vertices + 1
                temp_vertices(n_temp_vertices) = cell_vertex
              END IF
            END DO
          END DO

#endif

        END DO

        IF (n_temp_vertices .GT. 1) THEN
          CALL quicksort(temp_vertices(1:n_temp_vertices))
          CALL remove_duplicated_entries(temp_vertices(1:n_temp_vertices), &
                                       n_temp_vertices)
          ! remove vertices that are on level ilev and ilev - 1
          DO k = -1, 0
            CALL remove_entries_from_ref_list(temp_vertices(1:n_temp_vertices), &
              n_temp_vertices, flag2_v_list(ilev+k)%idx(1:n2_ilev_v(ilev+k)))
          END DO
        END IF

        IF (SIZE(temp_vertices_owner(:)) < n_temp_vertices) THEN
          DEALLOCATE(temp_vertices_owner)
          ALLOCATE(temp_vertices_owner(n_temp_vertices))
        END IF

        ALLOCATE(flag2_v_list(ilev+1)%idx(n_temp_vertices), &
          &      flag2_v_list(ilev+1)%owner(n_temp_vertices))
        flag2_v_list(ilev+1)%idx(:) = temp_vertices(1:n_temp_vertices)
        n2_ilev_v(ilev+1) = n_temp_vertices
        CALL compute_vertex_owner(flag2_v_list(ilev+1)%idx(:), &
          &                       flag2_v_list(ilev+1)%owner(:))
      END DO

      DEALLOCATE(temp_cells, temp_vertices, temp_edges, edge_cells, pack_mask)
    END SUBROUTINE compute_flag_lists

    SUBROUTINE compute_flag_lists_short(flag2_c_list, flag2_v_list, &
      &                                 flag2_e_list, n2_ilev_c, n2_ilev_v, &
      &                                 n2_ilev_e, n_boundary_rows &
    ! removed arguments because NEC has some problems with them...
    !  &                                 ,owned_edges, owned_verts &
                                        )
      INTEGER, INTENT(IN) :: n_boundary_rows
      TYPE(nb_flag_list_elem), INTENT(OUT) :: &
        flag2_c_list(0:2*n_boundary_rows), flag2_v_list(0:n_boundary_rows+1), &
        flag2_e_list(0:2*n_boundary_rows+1)
      INTEGER, INTENT(OUT) :: n2_ilev_c(0:2*n_boundary_rows), &
        &                     n2_ilev_v(0:n_boundary_rows+1), &
        &                     n2_ilev_e(0:2*n_boundary_rows+1)
      ! removed arguments because NEC has some problems with them...
      ! INTEGER, ALLOCATABLE, INTENT(OUT) :: owned_edges(:), owned_verts(:)

      INTEGER :: i
      LOGICAL :: position_error

      !---------------------------------------------------------------------
      ! flag_c_list(-1)%idx empty dummy list
      ! flag_c_list(0)%idx  all cells jg where cell_owner(jg) == my_proc
      ! flag_c_list(j)%idx j in 1..n_boundary_cells all cells bordering
      !                    on cells in flag_c_list(j - 1)
      ! flag_e_list(-1)%idx empty dummy list
      ! flag_e_list(j)%idx all edges where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the edge
      ! flag_v_list(-1)%idx empty dummy list
      ! flag_v_list(j)%idx all vertices where an adjacent cell is in
      !                    flag_c_list(j)%idx for minimal j,
      !                    i.e. no other cell in list flag_c_list(k)%idx
      !                    with k < j is adjacent to the vertex
      !------------------------------------------------------------------------
      ! The purpose of the second set of flags is to control moving the
      ! halo points to the end of the index vector for nearly all cells
      ! even if order_type_of_halos = 1
      !------------------------------------------------------------------------
      ! flag2_e_list(0)%idx all edges which are owned by the local process
      ! flag2_e_list(1)%idx all edges adjacent to cells in flag_c_list(0)
      !                     but owned by other tasks
      ! flag2_e_list(2*i)%idx where i > 0, edges in idx are ONLY adjacent
      !                       to cells in flag_c_list(i)%idx
      ! flag2_e_list(2*i+1)%idx where i > 0, edges adjacent to one cell in
      !                         flag_c_list(i) and either
      !                            one cell in flag_c_list(i+1)
      !                         or edge is an outer edge,
      !                            i.e. adjacent to only one cell
      ! flag2_c_list(0)%idx == flag_c_list(0)%idx
      ! flag2_c_list(2*j-1)%idx where j > 0, contains all cells from
      !                         flag_c_list(j) which are
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_c_list(2*j)%idx   where j > 0, contains all cells from
      !                         flag_c_list(j) which are NOT
      !                         adjacent to cells in flag_c_list(j-1)
      ! flag2_v_list(0)%idx contains all vertices from flag_v_list(0)%idx
      !                         owned by the current task
      ! flag2_v_list(1)%idx contains all vertices from flag_v_list(0)%idx
      !                     NOT owned by the current task
      ! flag2_v_list(j)%idx for j > 1 == flag_v_list(j-1)%idx
      !------------------------------------------------------------------------

      !--------------------------------------------------------------------------
      ! compute flag2 arrays
      !--------------------------------------------------------------------------

      position_error = .FALSE.
      DO i = 1, wrk_p_patch_pre%n_patch_cells_g
        position_error = (owned_cells(i) /= i) .OR. position_error
      END DO
      IF (position_error) &
        CALL finish("compute_flag_lists_short", &
          &         "owned_cells content is inconsistent")

      ! set n2_ilev_cve arrays
      n2_ilev_c(:) = 0
      n2_ilev_v(:) = 0
      n2_ilev_e(:) = 0
      n2_ilev_c(0) = wrk_p_patch_pre%n_patch_cells_g
      n2_ilev_v(0) = wrk_p_patch_pre%n_patch_verts_g
      n2_ilev_e(0) = wrk_p_patch_pre%n_patch_edges_g

      ! allocate and set flag2_cve_list arrays
      ALLOCATE(flag2_c_list(0)%idx(wrk_p_patch_pre%n_patch_cells_g), &
        &      flag2_c_list(0)%owner(wrk_p_patch_pre%n_patch_cells_g), &
        &      flag2_v_list(0)%idx(wrk_p_patch_pre%n_patch_verts_g), &
        &      flag2_v_list(0)%owner(wrk_p_patch_pre%n_patch_verts_g), &
        &      flag2_v_list(1)%idx(0), &
        &      flag2_v_list(1)%owner(0), &
        &      flag2_e_list(0)%idx(wrk_p_patch_pre%n_patch_edges_g), &
        &      flag2_e_list(0)%owner(wrk_p_patch_pre%n_patch_edges_g), &
        &      flag2_e_list(1)%idx(0), &
        &      flag2_e_list(1)%owner(0))

      DO i = 1, wrk_p_patch_pre%n_patch_cells_g
        flag2_c_list(0)%idx(i) = i
      END DO
      flag2_c_list(0)%owner(:) = my_proc

      DO i = 1, wrk_p_patch_pre%n_patch_verts_g
        flag2_v_list(0)%idx(i) = i
      END DO
      flag2_v_list(0)%owner(:) = my_proc

      DO i = 1, wrk_p_patch_pre%n_patch_edges_g
        flag2_e_list(0)%idx(i) = i
      END DO
      flag2_e_list(0)%owner(:) = my_proc

      DO i = 1, n_boundary_rows

        ALLOCATE(flag2_c_list(2*i-1)%idx(0), flag2_c_list(2*i)%idx(0), &
          &      flag2_c_list(2*i-1)%owner(0), flag2_c_list(2*i)%owner(0), &
          &      flag2_v_list(i+1)%idx(0), flag2_v_list(i+1)%owner(0), &
          &      flag2_e_list(2*i)%idx(0), flag2_e_list(2*i+1)%idx(0), &
          &      flag2_e_list(2*i)%owner(0), flag2_e_list(2*i+1)%owner(0))
      END DO

      ! allocate and set owner arrays
      ALLOCATE(owned_verts(wrk_p_patch_pre%n_patch_verts_g), &
        &      owned_edges(wrk_p_patch_pre%n_patch_edges_g))
      owned_verts(:) = flag2_v_list(0)%idx(:)
      owned_edges(:) = flag2_e_list(0)%idx(:)
    END SUBROUTINE compute_flag_lists_short

    SUBROUTINE compute_vertex_owner(vertices, owner)
      INTEGER, INTENT(IN) :: vertices(:)
      INTEGER, INTENT(OUT) :: owner(:)

      INTEGER :: i, jv, j, temp_num_owners, &
        &        a_iown(2), a_mod_iown(2)
      LOGICAL :: swap
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      INTEGER, ALLOCATABLE :: temp_num_edges(:), t_cells(:,:), t_cell_owner(:,:)
#else
      INTEGER :: t_cells(wrk_p_patch_pre%verts%max_connectivity), &
        &        t_cell_owner(wrk_p_patch_pre%verts%max_connectivity), &
        &        temp_num_edges
#endif
      ! compute ownership of vertices
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED

      ALLOCATE(temp_num_edges(SIZE(vertices)), &
        &      t_cells(SIZE(vertices),wrk_p_patch_pre%verts%max_connectivity), &
        &      t_cell_owner(SIZE(vertices),wrk_p_patch_pre%verts%max_connectivity))
#ifdef __BLOCK_GET__
      call dist_mult_init_blk_comm(1, SIZE(vertices), 1, wrk_p_patch_pre%verts%dist)
#endif
      DO i = 1, SIZE(vertices)
        jv = vertices(i)
#ifndef __BLOCK_GET__
        CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_num_edges, &
             (/jv/), temp_num_edges(i))
#else
        CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%verts%dist, v_num_edges, &
             (/jv/), temp_num_edges(i), proc_lst(1,i,1), ind_lst(1,i,1))
#endif
      END DO
#ifdef __BLOCK_GET__
      call dist_mult_do_blk_comm(wrk_p_patch_pre%verts%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%verts%dist)
#ifdef __BLOCK_GET__
      DO i = 1, SIZE(vertices)
        IF (proc_lst(1,i,1) /= -1) THEN
          temp_num_edges(i) = recv_array(proc_lst(1,i,1))%lst(ind_lst(1,i,1))
        END IF
      END DO

      call dist_mult_end_blk_comm(wrk_p_patch_pre%verts%dist)
      call dist_mult_init_blk_comm(1, SIZE(vertices), maxval(temp_num_edges(1:SIZE(vertices))), wrk_p_patch_pre%verts%dist)
#endif
      DO i = 1, SIZE(vertices)
        jv = vertices(i)
        DO j = 1, temp_num_edges(i)
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_cell, &
               (/jv, j/), t_cells(i,j))
#else
          CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%verts%dist, v_cell, &
               (/jv, j/), t_cells(i,j), proc_lst(1,i,j), ind_lst(1,i,j))
#endif
        END DO
      END DO
#ifdef __BLOCK_GET__
      call dist_mult_do_blk_comm(wrk_p_patch_pre%verts%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%verts%dist)
#ifdef __BLOCK_GET__
      DO i = 1, SIZE(vertices)
        DO j = 1, temp_num_edges(i)
          IF (proc_lst(1,i,j) /= -1) THEN
            t_cells(i,j) = recv_array(proc_lst(1,i,j))%lst(ind_lst(1,i,j))
          END IF
        END DO
      END DO

      call dist_mult_end_blk_comm(wrk_p_patch_pre%verts%dist)
      call dist_mult_init_blk_comm(1, SIZE(vertices), maxval(temp_num_edges(1:SIZE(vertices))), dist_cell_owner)
#endif
      DO i = 1, SIZE(vertices)
        CALL quicksort(t_cells(i,1:temp_num_edges(i)))
        DO j = 1, temp_num_edges(i)
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(dist_cell_owner, 1, (/t_cells(i,j)/), &
            &                      t_cell_owner(i,j))
#else
          CALL dist_mult_array_get_i4_blk(dist_cell_owner, 1, (/t_cells(i,j)/), &
            &                      t_cell_owner(i,j), proc_lst(1,i,j), ind_lst(1,i,j))
#endif
        END DO
      END DO
#ifdef __BLOCK_GET__
      call dist_mult_do_blk_comm(dist_cell_owner)
#endif
      CALL dist_mult_array_rma_sync(dist_cell_owner)

#ifdef __BLOCK_GET__
      DO i = 1, SIZE(vertices)
        DO j = 1, temp_num_edges(i)
          IF (proc_lst(1,i,j) /= -1) THEN
            t_cell_owner(i,j) = recv_array(proc_lst(1,i,j))%lst(ind_lst(1,i,j))
          END IF
        END DO
      END DO

      call dist_mult_end_blk_comm(dist_cell_owner)
#endif
      DO i = 1, SIZE(vertices)
        temp_num_owners = COUNT(t_cell_owner(i,1:temp_num_edges(i)) >= 0)
        t_cell_owner(i,1:temp_num_owners) = PACK(t_cell_owner(i,1:temp_num_edges(i)), &
          &                         t_cell_owner(i,1:temp_num_edges(i)) >= 0)
        a_iown(1) = t_cell_owner(i,1)
        a_mod_iown(1) = MOD(a_iown(1), 2)
        DO j = 2, temp_num_owners
          a_iown(2) = t_cell_owner(i,j)
          a_mod_iown(2) = MOD(a_iown(2), 2)
          ! 50% chance of acquiring a vertex
          swap = a_iown(1) > a_iown(2) &
               .NEQV. a_mod_iown(1) == a_mod_iown(2)
          a_iown(1) = MERGE(a_iown(2), a_iown(1), swap)
          a_mod_iown(1) = MERGE(a_mod_iown(2), a_mod_iown(1), swap)
        END DO
        owner(i) = a_iown(1)
      END DO
      DEALLOCATE(temp_num_edges, t_cells, t_cell_owner)

#else

      DO i = 1, SIZE(vertices)
        jv = vertices(i)
        CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_num_edges, &
             (/jv/), temp_num_edges)
        DO j = 1, temp_num_edges
          CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_cell, &
               (/jv, j/), t_cells(j))
        END DO
        IF (temp_num_edges .GT. 1) THEN
          CALL quicksort(t_cells(1:temp_num_edges))
        ENDIF
        DO j = 1, temp_num_edges
          CALL dist_mult_array_get(dist_cell_owner, 1, (/t_cells(j)/), &
            &                      t_cell_owner(j))
        END DO
        temp_num_owners = COUNT(t_cell_owner(1:temp_num_edges) >= 0)
        t_cell_owner(1:temp_num_owners) = PACK(t_cell_owner(1:temp_num_edges), &
             t_cell_owner(1:temp_num_edges) >= 0)
        a_iown(1) = t_cell_owner(1)
        a_mod_iown(1) = MOD(a_iown(1), 2)
        DO j = 2, temp_num_owners
          a_iown(2) = t_cell_owner(j)
          a_mod_iown(2) = MOD(a_iown(2), 2)
          ! 50% chance of acquiring a vertex
          swap = a_iown(1) > a_iown(2) &
               .NEQV. a_mod_iown(1) == a_mod_iown(2)
          a_iown(1) = MERGE(a_iown(2), a_iown(1), swap)
          a_mod_iown(1) = MERGE(a_mod_iown(2), a_mod_iown(1), swap)
        END DO
        owner(i) = a_iown(1)
      END DO

#endif

    END SUBROUTINE compute_vertex_owner

    SUBROUTINE compute_edge_owner(edges, owner)
      INTEGER, INTENT(IN) :: edges(:)
      INTEGER, INTENT(OUT) :: owner(:)

      INTEGER :: i, je, a_mod_iown(2)

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      INTEGER, ALLOCATABLE :: a_idx(:,:), a_iown(:,:)
#else
      INTEGER :: a_idx(2), a_iown(2)
#endif

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      ALLOCATE(a_idx(SIZE(edges), 2), a_iown(SIZE(edges), 2))
#ifdef __BLOCK_GET__
      call dist_mult_init_blk_comm(2, SIZE(edges),1,wrk_p_patch_pre%edges%dist)
#endif
      DO i = 1, SIZE(edges)

        je = edges(i)
#ifndef __BLOCK_GET__
        CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
             (/je,1/), a_idx(i, 1))
        CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
             (/je,2/),a_idx(i, 2))
#else
        CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_cell, &
             (/je,1/), a_idx(i, 1), proc_lst(1,i,1), ind_lst(1,i,1))
        CALL dist_mult_array_get_i4_blk(wrk_p_patch_pre%edges%dist, e_cell, &
             (/je,2/), a_idx(i, 2), proc_lst(2,i,1), ind_lst(2,i,1))
#endif
      END DO
#ifdef __BLOCK_GET__
      CALL dist_mult_do_blk_comm(wrk_p_patch_pre%edges%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%edges%dist)
#ifdef __BLOCK_GET__
      DO i = 1, SIZE(edges)
        IF (proc_lst(1,i,1) /= -1) THEN
          a_idx(i, 1) = recv_array(proc_lst(1,i,1))%lst(ind_lst(1,i,1))
        END IF
        IF (proc_lst(2,i,1) /= -1) THEN
          a_idx(i, 2) = recv_array(proc_lst(2,i,1))%lst(ind_lst(2,i,1))
        END IF
      END DO
      CALL dist_mult_end_blk_comm(wrk_p_patch_pre%edges%dist)
#endif

#ifdef __BLOCK_GET__
      CALL dist_mult_init_blk_comm(2, SIZE(edges),1,dist_cell_owner)
#endif
      DO i = 1, SIZE(edges)
        ! outer boundary edges always belong to single adjacent cell owner
        a_idx(i,:) = MERGE(a_idx(i,:), MAXVAL(a_idx(i,:)), a_idx(i,:) > 0)
#ifndef __BLOCK_GET__
        CALL dist_mult_array_get(dist_cell_owner, 1, (/a_idx(i,1)/), a_iown(i,1))
        CALL dist_mult_array_get(dist_cell_owner, 1, (/a_idx(i,2)/), a_iown(i,2))
#else
        CALL dist_mult_array_get_i4_blk(dist_cell_owner, 1, (/a_idx(i,1)/), a_iown(i,1), &
             proc_lst(1,i,1), ind_lst(1,i,1))
        CALL dist_mult_array_get_i4_blk(dist_cell_owner, 1, (/a_idx(i,2)/), a_iown(i,2), &
             proc_lst(2,i,1), ind_lst(2,i,1))
#endif
      END DO
#ifdef __BLOCK_GET__
      CALL dist_mult_do_blk_comm(dist_cell_owner)
#endif
      CALL dist_mult_array_rma_sync(dist_cell_owner)
#ifdef __BLOCK_GET__
      DO i = 1, SIZE(edges)
        IF (proc_lst(1,i,1) /= -1) THEN
          a_iown(i,1) = recv_array(proc_lst(1,i,1))%lst(ind_lst(1,i,1))
        END IF
        IF (proc_lst(2,i,1) /= -1) THEN
          a_iown(i,2) = recv_array(proc_lst(2,i,1))%lst(ind_lst(2,i,1))
        END IF
      END DO
      CALL dist_mult_end_blk_comm(dist_cell_owner)
#endif
      DO i = 1, SIZE(edges)

        a_iown(i,:) = MERGE(a_iown(i,:), MAXVAL(a_iown(i,:)), a_iown(i,:) >= 0)
        a_mod_iown(1:2) = MOD(a_iown(i,1:2), 2)
        ! 50% chance of acquiring an edge
        owner(i) = MERGE(a_iown(i,1), a_iown(i,2), &
          &              (a_iown(i,1) > a_iown(i,2)) .EQV. &
          &              (a_mod_iown(1) == a_mod_iown(2)))
      END DO
      DEALLOCATE(a_idx, a_iown)

#else

      DO i = 1, SIZE(edges)

        je = edges(i)
        CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
             (/je,1/), a_idx(1))
        CALL dist_mult_array_get(wrk_p_patch_pre%edges%dist, e_cell, &
             (/je,2/), a_idx(2))
        ! outer boundary edges always belong to single adjacent cell owner
        a_idx(:) = MERGE(a_idx(:), MAXVAL(a_idx(:)), a_idx(:) > 0)
        CALL dist_mult_array_get(dist_cell_owner, 1, (/a_idx(1)/), a_iown(1))
        CALL dist_mult_array_get(dist_cell_owner, 1, (/a_idx(2)/), a_iown(2))
        a_iown(:) = MERGE(a_iown(:), MAXVAL(a_iown(:)), a_iown(:) >= 0)
        a_mod_iown(1:2) = MOD(a_iown(1:2), 2)
        ! 50% chance of acquiring an edge
        owner(i) = MERGE(a_iown(1), a_iown(2), &
          &              (a_iown(1) > a_iown(2)) .EQV. &
          &              (a_mod_iown(1) == a_mod_iown(2)))
      END DO

#endif

    END SUBROUTINE compute_edge_owner

    SUBROUTINE sort_flag_list_by_owner(flag_list)
      TYPE(nb_flag_list_elem), INTENT(INOUT) :: flag_list
      INTEGER :: i, prev_owner, prev_pos

      IF (SIZE(flag_list%owner) <= 1) RETURN

      CALL quicksort(flag_list%owner, flag_list%idx)

      prev_owner = flag_list%owner(1)
      prev_pos = 1
      DO i = 2, SIZE(flag_list%owner)
        IF (flag_list%owner(i) /= prev_owner) THEN
          CALL quicksort(flag_list%idx(prev_pos:i-1))
          prev_pos = i
          prev_owner = flag_list%owner(i)
        END IF
      END DO

      CALL quicksort(flag_list%idx(prev_pos:))
    END SUBROUTINE

  END SUBROUTINE divide_patch

  FUNCTION refine_cells(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag == 0 &
         .OR. ( flag >= 1 .AND. irl0 == 1 ) &
         .OR. ( flag == 1 .AND. irl0 == 2 )
  END FUNCTION refine_cells

  FUNCTION refine_edges(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag==0 &
         .OR. ((irl0==1 .OR. irl0==2) .AND. flag>=1) &
         .OR. ((irl0==3 .OR. irl0==4) .AND. (flag==1 .OR. flag==2)) &
         .OR. ((irl0==5 .OR. irl0==6 .OR. irl0<0) .AND. flag==1)
  END FUNCTION refine_edges

  FUNCTION refine_verts(flag, irl0) RESULT(p)
    INTEGER, INTENT(in) :: flag, irl0
    LOGICAL :: p

    p = flag == 0 &
         .OR.  (irl0 == 1                .AND. flag  > 0) &
         .OR. ((irl0 == 2 .OR. irl0 < 0) .AND. flag == 1)
  END FUNCTION refine_verts

  SUBROUTINE build_patch_start_end(n_patch_cve, n_patch_cve_g, &
       patch_id, nblks, npromz, max_cell_connectivity, &
       min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve, &
       flag2_list, n2_ilev, decomp_info, &
       start_index, start_block, end_index, end_block, &
       start_g, end_g, order_type_of_halos, l_cell_correction, refin_ctrl, &
       refin_ctrl_aidx, refinement_predicate)
    INTEGER, INTENT(in) :: n_patch_cve, n_patch_cve_g, &
         patch_id, nblks, npromz, max_cell_connectivity, &
         min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve
    INTEGER, INTENT(in) :: order_type_of_halos, refin_ctrl_aidx
    LOGICAL, INTENT(in) :: l_cell_correction
    TYPE(dist_mult_array) :: refin_ctrl
    INTEGER, INTENT(in) :: n2_ilev(0:)
    TYPE(nb_flag_list_elem), INTENT(in) :: flag2_list(0:)
    TYPE(t_grid_domain_decomp_info), INTENT(inout) :: decomp_info
    INTEGER, DIMENSION(min_rlcve:), INTENT(inout) :: &
         start_index, start_block, end_index, end_block
    INTEGER, DIMENSION(min_rlcve:), INTENT(in) :: start_g, end_g
    INTERFACE
      FUNCTION refinement_predicate(flag, irl0) RESULT(p)
        INTEGER, INTENT(in) :: flag, irl0
        LOGICAL :: p
      END FUNCTION refinement_predicate
    END INTERFACE

    CHARACTER(LEN=*), PARAMETER :: routine = modname//'::build_patch_start_end'

    INTEGER :: i, ilev, ilev1, ilev_st, irlev, j, jb, jf, jl, &
         exec_sequence, ref_flag, k, k_start, n_inner, temp
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    INTEGER, ALLOCATABLE, DIMENSION(:) :: irl0
#else
    INTEGER :: irl0
#endif
    INTEGER, ALLOCATABLE :: temp_glb_index(:), permutation(:), temp_ilev(:), &
      &                     temp_owner(:), temp_copy(:)

    ! if all cells/vertices/edges have flag == 0
    IF ((n2_ilev(0) == n_patch_cve) .AND. (n2_ilev(0) == n_patch_cve_g)) THEN

      CALL build_patch_start_end_short(n_patch_cve, n_patch_cve_g, &
        patch_id, nblks, npromz, max_cell_connectivity, &
        min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve, &
        flag2_list, n2_ilev, decomp_info, &
        start_index, start_block, end_index, end_block, start_g, end_g, &
        order_type_of_halos, l_cell_correction)
      RETURN
    END IF

    ALLOCATE(temp_glb_index(n_patch_cve), temp_owner(n_patch_cve))

    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      decomp_info%glb_index(1:n2_ilev(0)) = flag2_list(0)%idx(1:n2_ilev(0))
      decomp_info%owner_local(1:n2_ilev(0)) = flag2_list(0)%owner(1:n2_ilev(0))
      n_inner = n2_ilev(0)

    CASE (1)
      ALLOCATE(temp_ilev(n_patch_cve), permutation(n_patch_cve), &
           temp_copy(n_patch_cve))
      k = 1
      DO ilev = 0, max_ilev
        temp_glb_index(k:k+n2_ilev(ilev)-1) = flag2_list(ilev)%idx(1:n2_ilev(ilev))
        temp_copy(k:k+n2_ilev(ilev)-1) = flag2_list(ilev)%owner(1:n2_ilev(ilev))
        temp_ilev(k:k+n2_ilev(ilev)-1) = ilev
        k = k + n2_ilev(ilev)
      END DO
      DO k = 1, n_patch_cve
        permutation(k) = k
      END DO
      CALL quicksort(temp_glb_index, permutation)
      DO k = 1, n_patch_cve
        temp_owner(k) = temp_copy(permutation(k))
      END DO
      temp_copy = temp_ilev
      DO k = 1, n_patch_cve
        temp_ilev(k) = temp_copy(permutation(k))
      END DO
      DEALLOCATE(permutation, temp_copy)
      k = 1
      jf = 1
      j = 0
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      ALLOCATE(irl0(n_patch_cve))
      DO WHILE (j < n_patch_cve_g .AND. jf <= n_patch_cve)
        j = temp_glb_index(jf)
        CALL dist_mult_array_get(refin_ctrl, refin_ctrl_aidx, (/j/), irl0(jf))
        jf = jf + 1
      END DO
      CALL dist_mult_array_rma_sync(refin_ctrl)
      DO i = 1, jf - 1
        IF (refinement_predicate(temp_ilev(i), irl0(i))) THEN
          decomp_info%glb_index(k) = temp_glb_index(i)
          decomp_info%owner_local(k) = temp_owner(i)
          k = k + 1
        END IF
      END DO
      DEALLOCATE(irl0)
#else
      DO WHILE (j < n_patch_cve_g .AND. jf <= n_patch_cve)
        j = temp_glb_index(jf)
        CALL dist_mult_array_get(refin_ctrl, refin_ctrl_aidx, (/j/), irl0)
        IF (refinement_predicate(temp_ilev(jf), irl0)) THEN
          decomp_info%glb_index(k) = j
          decomp_info%owner_local(k) = temp_owner(jf)
          k = k + 1
        END IF
        jf = jf + 1
      END DO
#endif
      n_inner = k - 1
      DEALLOCATE(temp_ilev, temp_glb_index, temp_owner)

    CASE default
      CALL finish("", "Uknown order_type_of_halos")
    END SELECT

    ALLOCATE(permutation(n_inner))
    DO k = 1, n_inner
      permutation(k) = k
    END DO
    CALL set_inner_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(1:n_inner), &
      &                      permutation)
    DEALLOCATE(permutation)
    ! Set start_index/block ... end_index/block for cells/verts/edges.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the inner global indices.

    ! Note: the segments between min_rlcell and min_rlcell_int-1 are reserved for
    ! halo cells; they are set below
    DO i=min_rlcve_int,max_rlcve
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, start_g(i), temp, +1)
      start_index(i) = idx_no(temp)
      start_block(i) = blk_no(temp)
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, end_g(i), temp, -1)
      end_index(i) = idx_no(temp)
      end_block(i) = blk_no(temp)
    ENDDO

    ! Preset remaining index sections with dummy values
    start_index(min_rlcve:min_rlcve_int-1) = -9999
    start_block(min_rlcve:min_rlcve_int-1) = -9999
    end_index(min_rlcve:min_rlcve_int-1)   = -9999
    end_block(min_rlcve:min_rlcve_int-1)   = -9999

    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      ! processing local parent grid
      !
      ! Gather all halo cells/edges/verts at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rlcve_int-1
      irlev = min_rlcve_int-1
      ref_flag = MERGE(0, 1, l_cell_correction .OR. order_type_of_halos == 2)
      ! write(0,*) "ref_flag=", ref_flag
      j = 1
      DO ilev = ref_flag + 1, max_ilev
        temp_glb_index(j:j+n2_ilev(ilev)-1) &
             = flag2_list(ilev)%idx(1:n2_ilev(ilev))
        temp_owner(j:j+n2_ilev(ilev)-1) &
             = flag2_list(ilev)%owner(1:n2_ilev(ilev))
        j = j + n2_ilev(ilev)
      END DO
      CALL quicksort(temp_glb_index(1:j-1), temp_owner(1:j-1))
      k_start = SUM(n2_ilev(0:ref_flag))
      IF (start_index(irlev)==-9999 .AND. k_start + 1 <= n_patch_cve) THEN
        start_index(irlev) = idx_no(k_start + 1) ! line index
        start_block(irlev) = blk_no(k_start + 1) ! block index
      ENDIF

      DO k = k_start + 1, n_patch_cve
        decomp_info%glb_index(k) = temp_glb_index(k - k_start)
        decomp_info%owner_local(k) = temp_owner(k - k_start)
      END DO
      DEALLOCATE(temp_glb_index, temp_owner)

      ! Set end index when the last point is processed
      end_index(min_rlcve:irlev) = idx_no(n_patch_cve)
      end_block(min_rlcve:irlev) = blk_no(n_patch_cve)
      start_index(min_rlcve:irlev-1) = start_index(irlev)
      start_block(min_rlcve:irlev-1) = start_block(irlev)

    CASE(1)
      k = n_inner
      ! Gather all halo cells/verts/edges except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2 value and will be placed at the
      ! index levels between min_rlcve_int-1 and min_rlcve
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here
      DO ilev = 1, max_ilev
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
        ALLOCATE(irl0(n2_ilev(ilev)))
#ifdef __BLOCK_GET__
        CALL dist_mult_init_blk_comm(1, n2_ilev(ilev),1,refin_ctrl)
#endif
        DO j = 1, n2_ilev(ilev)
          jf = flag2_list(ilev)%idx(j)
#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(refin_ctrl, refin_ctrl_aidx, (/jf/), irl0(j))
#else
          CALL dist_mult_array_get_i4_blk(refin_ctrl, refin_ctrl_aidx, (/jf/), irl0(j), &
               proc_lst(1,j,1), ind_lst(1,j,1))
#endif
        END DO
#ifdef __BLOCK_GET__
        CALL dist_mult_do_blk_comm(refin_ctrl)
#endif
        CALL dist_mult_array_rma_sync(refin_ctrl)
#ifdef __BLOCK_GET__
        DO j = 1, n2_ilev(ilev)
          IF (proc_lst(1,j,1) /= -1) THEN
            irl0(j) = recv_array(proc_lst(1,j,1))%lst(ind_lst(1,j,1))
          END IF
        END DO
        CALL dist_mult_end_blk_comm(refin_ctrl)
#endif
        DO j = 1, n2_ilev(ilev)
          jf = flag2_list(ilev)%idx(j)
          IF (.NOT. refinement_predicate(ilev, irl0(j))) THEN
            k = k + 1
            decomp_info%glb_index(k) = jf
            decomp_info%owner_local(k) = flag2_list(ilev)%owner(j)
            jb = blk_no(k) ! block index
            jl = idx_no(k) ! line index
            ! This ensures that just the first point found at this irlev is saved as start point
            IF (start_index(irlev)==-9999) THEN
              start_index(irlev) = jl
              start_block(irlev) = jb
            ENDIF
            end_index(irlev) = jl
            end_block(irlev) = jb
          END IF
        END DO
        DEALLOCATE(irl0)
#else
        DO j = 1, n2_ilev(ilev)
          jf = flag2_list(ilev)%idx(j)
          CALL dist_mult_array_get(refin_ctrl, refin_ctrl_aidx, (/jf/), irl0)
          IF (.NOT. refinement_predicate(ilev, irl0)) THEN
            k = k + 1
            decomp_info%glb_index(k) = jf
            decomp_info%owner_local(k) = flag2_list(ilev)%owner(j)
            jb = blk_no(k) ! block index
            jl = idx_no(k) ! line index
            ! This ensures that just the first point found at this irlev is saved as start point
            IF (start_index(irlev)==-9999) THEN
              start_index(irlev) = jl
              start_block(irlev) = jb
            ENDIF
            end_index(irlev) = jl
            end_block(irlev) = jb
          END IF
        END DO
#endif

        ! Just in case that no grid point is found (may happen for ilev=1)
        IF (start_index(irlev)==-9999) THEN
          start_index(irlev) = end_index(irlev+1)+1
          start_block(irlev) = end_block(irlev+1)
          end_index(irlev)   = end_index(irlev+1)
          end_block(irlev)   = end_block(irlev+1)
        ENDIF

      ENDDO

      ! Fill start and end indices for remaining index sections
      IF (patch_id == 0) THEN
        ilev1 = min_rlcve_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rlcve, min_rlcve_int - max_ilev)
        ilev_st = max_ilev + 1
      ENDIF
      IF (l_cell_correction .AND. max_cell_connectivity==6) THEN ! for hexagons, there are no even-order halo cells
        DO ilev = 2, max_hw_cve, 2
          irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
          start_index(irlev) = end_index(irlev+1) + 1
          start_block(irlev) = end_block(irlev+1)
          end_index(irlev)   = end_index(irlev+1)
          end_block(irlev)   = end_block(irlev+1)
        ENDDO
      ENDIF
      DO ilev = ilev_st,max_hw_cve
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
        start_index(irlev) = end_index(ilev1) + 1
        start_block(irlev) = end_block(ilev1)
        end_index(irlev)   = end_index(ilev1)
        end_block(irlev)   = end_block(ilev1)
      ENDDO
    END SELECT

    ! exec_sequence only serves the purpose to retain the execution
    ! sequence for these three correction as it was in the previous
    ! code version, it might be irrelevant.
    ! see redmine issue #3924
    DO exec_sequence = 0, 1
      IF (.NOT. l_cell_correction .AND. exec_sequence == 0 &
           .OR. l_cell_correction .AND. exec_sequence == 1) THEN
        ! If a PE owns only nest boundary points, it may happen that
        ! one or more index sections of halo cells/verts/edges are
        ! empty. These are filled here
        IF (start_block(0)     > end_block(min_rlcve_int) &
             .OR. start_block(0) == end_block(min_rlcve_int) &
             .AND. start_index(0) > end_index(min_rlcve_int) ) THEN
          DO i = min_rlcve_int-1, min_rlcve, -1
            IF (start_index(i) == -9999 .OR. &
                 start_block(i) == -9999 .OR. &
                 end_index(i)   == -9999 .OR. &
                 end_block(i)   == -9999 ) THEN
              start_index(i) = start_index(i+1)
              start_block(i) = start_block(i+1)
              end_index(i)   = end_index(i+1)
              end_block(i)   = end_block(i+1)
            ENDIF
          ENDDO
        ENDIF

        ! Finally, fill start indices of halo rows with meaningful
        ! values for empty patches (which occur when processor
        ! splitting is applied)
        IF (nblks <= 0 .OR. npromz <= 0) THEN
          DO i=min_rlcve,min_rlcve_int-1
            start_index(i) = start_index(min_rlcve_int)
            start_block(i) = start_block(min_rlcve_int)
            end_index(i)   = end_index(min_rlcve_int)
            end_block(i)   = end_block(min_rlcve_int)
          ENDDO
        ENDIF
      END IF
      ! special treatment for sequential runs with trivial decomposition
      IF (exec_sequence == 0 .AND. SUM(n2_ilev(1:)) == 0) THEN
        end_index(min_rlcve:min_rlcve_int-1) &
             = end_index(min_rlcve_int)
        end_block(min_rlcve:min_rlcve_int-1) &
             = end_block(min_rlcve_int)
        IF (end_index(min_rlcve_int) == nproma) THEN
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int) + 1
          start_index(min_rlcve:min_rlcve_int-1) = 1
        ELSE
          start_index(min_rlcve:min_rlcve_int-1) = &
               &    end_index(min_rlcve_int) + 1
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int)
        END IF
      END IF

    END DO
    ALLOCATE(permutation(n_inner+1:n_patch_cve))
    DO k = n_inner+1, n_patch_cve
      permutation(k) = k
    END DO
    CALL set_outer_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(n_inner+1:), &
      &                      permutation)
  END SUBROUTINE build_patch_start_end

  SUBROUTINE build_patch_start_end_short(n_patch_cve, n_patch_cve_g, &
       patch_id, nblks, npromz, max_cell_connectivity, &
       min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve, &
       flag2_list, n2_ilev, decomp_info, &
       start_index, start_block, end_index, end_block, start_g, end_g, &
       order_type_of_halos, l_cell_correction)
    INTEGER, INTENT(in) :: n_patch_cve, n_patch_cve_g, &
         patch_id, nblks, npromz, max_cell_connectivity, &
         min_rlcve, min_rlcve_int, max_rlcve, max_ilev, max_hw_cve
    INTEGER, INTENT(in) :: order_type_of_halos
    LOGICAL, INTENT(in) :: l_cell_correction
    INTEGER, INTENT(in) :: n2_ilev(0:)
    TYPE(nb_flag_list_elem), INTENT(in) :: flag2_list(0:)
    TYPE(t_grid_domain_decomp_info), INTENT(inout) :: decomp_info
    INTEGER, DIMENSION(min_rlcve:), INTENT(inout) :: &
         start_index, start_block, end_index, end_block
    INTEGER, DIMENSION(min_rlcve:), INTENT(in) :: start_g, end_g

    INTEGER :: i, ilev, ilev1, ilev_st, irlev, exec_sequence, temp

    ! if all cells/vertices/edges have flag == 0
    IF ((n2_ilev(0) /= n_patch_cve) .OR. (n2_ilev(0) /= n_patch_cve_g)) &
      CALL finish("build_patch_start_end_short", &
          &       "arguments do not fit for this routine")

    IF (ANY(order_type_of_halos == (/0,1,2/))) THEN
      decomp_info%glb_index(:) = flag2_list(0)%idx(:)
      decomp_info%owner_local(:) = flag2_list(0)%owner(:)
    ELSE
      CALL finish("build_patch_start_end_short", "Uknown order_type_of_halos")
    END IF

    CALL set_inner_glb_index(decomp_info%glb2loc_index, &
      &                      decomp_info%glb_index(:), &
      &                      (/(i, i = 1, n_patch_cve)/))

    ! Set start_index/block ... end_index/block for cells/verts/edges.
    ! This must be done here since it depends on the special (monotonic)
    ! setting of the inner global indices.

    ! Note: the segments between min_rlcell and min_rlcell_int-1 are reserved for
    ! halo cells; they are set below
    DO i=min_rlcve_int,max_rlcve
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, start_g(i), temp, +1)
      start_index(i) = idx_no(temp)
      start_block(i) = blk_no(temp)
!CDIR IEXPAND
      CALL get_local_idx(decomp_info, end_g(i), temp, -1)
      end_index(i) = idx_no(temp)
      end_block(i) = blk_no(temp)
    ENDDO

    ! Preset remaining index sections with dummy values
    start_index(min_rlcve:min_rlcve_int-1) = -9999
    start_block(min_rlcve:min_rlcve_int-1) = -9999
    end_index(min_rlcve:min_rlcve_int-1)   = -9999
    end_block(min_rlcve:min_rlcve_int-1)   = -9999

    SELECT CASE(order_type_of_halos)
    CASE (0,2)
      ! processing local parent grid
      !
      ! Gather all halo cells/edges/verts at the end of the patch.
      ! They do not undergo further ordering and will be placed at index level min_rlcve_int-1
      irlev = min_rlcve_int-1

      ! Set end index when the last point is processed
      start_index(min_rlcve:irlev-1) = start_index(irlev)
      start_block(min_rlcve:irlev-1) = start_block(irlev)
      end_index(min_rlcve:irlev) = idx_no(n_patch_cve)
      end_block(min_rlcve:irlev) = blk_no(n_patch_cve)

    CASE(1)
      ! Gather all halo cells/verts/edges except those lying in the lateral boundary interpolation zone
      ! at the end. They are sorted by the flag2 value and will be placed at the
      ! index levels between min_rlcve_int-1 and min_rlcve
      ! Note: this index range is empty on exit of prepare_gridref; therefore, get_local_index
      ! cannot be used here

      irlev = MAX(min_rlcve, min_rlcve_int - 1)  ! index section into which the halo points are put

      ! Just in case that no grid point is found (may happen for ilev=1)
      IF (.NOT. l_cell_correction) THEN
        start_index(irlev) = end_index(irlev+1)+1
        start_block(irlev) = end_block(irlev+1)
        end_index(irlev)   = end_index(irlev+1)
        end_block(irlev)   = end_block(irlev+1)
      ENDIF

      ! Fill start and end indices for remaining index sections
      IF (patch_id == 0) THEN
        ilev1 = min_rlcve_int
        ilev_st = 1
      ELSE
        ilev1 = MAX(min_rlcve, min_rlcve_int - max_ilev)
        ilev_st = max_ilev + 1
      ENDIF
      IF (l_cell_correction .AND. max_cell_connectivity==6) THEN ! for hexagons, there are no even-order halo cells
        DO ilev = 2, max_hw_cve, 2
          irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
          start_index(irlev) = end_index(irlev+1) + 1
          start_block(irlev) = end_block(irlev+1)
          end_index(irlev)   = end_index(irlev+1)
          end_block(irlev)   = end_block(irlev+1)
        ENDDO
      ENDIF
      DO ilev = ilev_st,max_hw_cve
        irlev = MAX(min_rlcve, min_rlcve_int - ilev)  ! index section into which the halo points are put
        start_index(irlev) = end_index(ilev1) + 1
        start_block(irlev) = end_block(ilev1)
        end_index(irlev)   = end_index(ilev1)
        end_block(irlev)   = end_block(ilev1)
      ENDDO
    END SELECT

    ! exec_sequence only serves the purpose to retain the execution
    ! sequence for these three correction as it was in the previous
    ! code version, it might be irrelevant.
    ! see redmine issue #3924
    DO exec_sequence = 0, 1
      IF (.NOT. l_cell_correction .AND. exec_sequence == 0 &
           .OR. l_cell_correction .AND. exec_sequence == 1) THEN
        ! If a PE owns only nest boundary points, it may happen that
        ! one or more index sections of halo cells/verts/edges are
        ! empty. These are filled here
        IF (start_block(0)     > end_block(min_rlcve_int) &
             .OR. start_block(0) == end_block(min_rlcve_int) &
             .AND. start_index(0) > end_index(min_rlcve_int) ) THEN
          DO i = min_rlcve_int-1, min_rlcve, -1
            IF (start_index(i) == -9999 .OR. &
                 start_block(i) == -9999 .OR. &
                 end_index(i)   == -9999 .OR. &
                 end_block(i)   == -9999 ) THEN
              start_index(i) = start_index(i+1)
              start_block(i) = start_block(i+1)
              end_index(i)   = end_index(i+1)
              end_block(i)   = end_block(i+1)
            ENDIF
          ENDDO
        ENDIF

        ! Finally, fill start indices of halo rows with meaningful
        ! values for empty patches (which occur when processor
        ! splitting is applied)
        IF (nblks <= 0 .OR. npromz <= 0) THEN
          DO i=min_rlcve,min_rlcve_int-1
            start_index(i) = start_index(min_rlcve_int)
            start_block(i) = start_block(min_rlcve_int)
            end_index(i)   = end_index(min_rlcve_int)
            end_block(i)   = end_block(min_rlcve_int)
          ENDDO
        ENDIF
      END IF
      ! special treatment for sequential runs with trivial decomposition
      IF (exec_sequence == 0 .AND. SUM(n2_ilev(1:)) == 0) THEN
        end_index(min_rlcve:min_rlcve_int-1) &
             = end_index(min_rlcve_int)
        end_block(min_rlcve:min_rlcve_int-1) &
             = end_block(min_rlcve_int)
        IF (end_index(min_rlcve_int) == nproma) THEN
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int) + 1
          start_index(min_rlcve:min_rlcve_int-1) = 1
        ELSE
          start_index(min_rlcve:min_rlcve_int-1) = &
               &    end_index(min_rlcve_int) + 1
          start_block(min_rlcve:min_rlcve_int-1) = &
               &   end_block(min_rlcve_int)
        END IF
      END IF
    END DO

  END SUBROUTINE build_patch_start_end_short

  !-------------------------------------------------------------------------
  !! Calculates local indices l_i from global indices g_i
  !! using the mapping in decomp_info
  !!
  SUBROUTINE get_local_idx(decomp_info, g_i, l_i, opt_mode)

    !
    TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
    INTEGER, INTENT(in) :: g_i
    INTEGER, INTENT(in), OPTIONAL :: opt_mode

    INTEGER, INTENT(out) :: l_i

    INTEGER :: mode

    !-----------------------------------------------------------------------

    ! mode controls the behaviour if the global index is valid, but outside the local domain:
    ! mode > 0 : Map it to the next greater index (even if this is outside the local domain)
    ! mode < 0 : Map it to the next smaller index (even if this is 0)
    ! mode = 0 : Map it to a negative value (negative of global index)

    IF(PRESENT(opt_mode)) THEN
      mode = opt_mode
    ELSE
      mode = 0
    ENDIF

    IF (mode == 0) THEN
      l_i = get_local_index(decomp_info%glb2loc_index, g_i)
      IF (l_i < 0) l_i = -g_i
    ELSE IF (mode > 0) THEN
      l_i = get_valid_local_index(decomp_info%glb2loc_index, g_i, .TRUE.)
    ELSE
      l_i = get_valid_local_index(decomp_info%glb2loc_index, g_i)
    END IF

  END SUBROUTINE get_local_idx

  !-------------------------------------------------------------------------
  !! Calculates local line/block indices l_idx, l_blk
  !! from global line/block indices g_idx, g_blk
  !! using the mapping in decomp_info
  !!
  SUBROUTINE get_local_idx_blk(decomp_info, g_idx, g_blk, l_idx, l_blk, opt_mode)

    !
    TYPE(t_grid_domain_decomp_info), INTENT(in) :: decomp_info
    INTEGER, INTENT(in) :: g_idx, g_blk
    INTEGER, INTENT(in), OPTIONAL :: opt_mode

    INTEGER, INTENT(out) :: l_idx, l_blk

    INTEGER j_l

    CALL get_local_idx(decomp_info, idx_1d(g_idx, g_blk), j_l, opt_mode)

    l_idx = idx_no(j_l)
    l_blk = blk_no(j_l)

  END SUBROUTINE get_local_idx_blk
  !-------------------------------------------------------------------------
  !! Makes a area subdivision for a subset of wrk_p_patch.
  !!
  SUBROUTINE divide_subset_geometric(subset_flag, n_proc, wrk_p_patch_pre, &
                                     dist_cell_owner, lparent_level, &
                                     divide_for_radiation)

    INTEGER, ALLOCATABLE   :: subset_flag(:) ! if > 0 a cell belongs to the subset
    INTEGER, INTENT(in)    :: n_proc   ! Number of processors
    TYPE(t_pre_patch), INTENT(INOUT) :: wrk_p_patch_pre ! out is required for
                                                        ! the distributed arrays
    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner ! receives the owner PE for every cell
    ! (-1 for cells not in subset)
    LOGICAL, INTENT(IN)    :: lparent_level ! indicates if domain decomposition is executed for
                              ! the parent grid level (true) or for the target grid level (false)
    ! Private flag if patch should be divided for radiation calculation
    LOGICAL, INTENT(in) :: divide_for_radiation

    INTEGER :: j, nc, nc_g, ncs, nce, n_onb_points, jm(1)
    INTEGER :: count_physdom(0:max_phys_dom), count_total, id_physdom(max_phys_dom), &
               num_physdom, proc_count(max_phys_dom), proc_offset(max_phys_dom), checksum, &
               ncell_offset(0:max_phys_dom), ncell_offset_g(0:max_phys_dom)
    TYPE(t_cell_info), ALLOCATABLE :: cell_desc(:)
    REAL(wp) :: corr_ratio(max_phys_dom)
    LOGICAL  :: lsplit_merged_domains, locean
    INTEGER :: my_cell_start, my_cell_end
    INTEGER, POINTER :: owner_chunk_ptr(:)
    INTEGER :: ierror
    CHARACTER(*), PARAMETER :: routine = 'divide_subset_geometric'

    !-----------------------------------------------------------------------

    IF(p_pe_work==0) THEN
      IF(divide_for_radiation) THEN
        WRITE(0,*) 'divide_patch: Using geometric area subdivision for radiation'
      ELSE
        WRITE(0,*) 'divide_patch: Using geometric area subdivision (normal)'
      ENDIF
    ENDIF

    locean = my_process_is_oceanic().or.my_process_is_waves()

    lsplit_merged_domains = .FALSE.

    my_cell_start = LBOUND(subset_flag, 1)
    my_cell_end = UBOUND(subset_flag, 1)
    IF (p_n_work > 1) THEN
#ifndef NOMPI
      CALL init_divide_cells_by_location_mpi
#else
      CALL finish("number of processors > 1 but MPI unavailable!")
#endif
    END IF

    ! Check if domain merging has been applied; for large PE numbers, calculating the DD
    ! for each physical domain separately tends to yield a more balanced distribution
    ! of the halo points
    IF (ldiv_phys_dom .AND. .NOT. divide_for_radiation) THEN
      count_physdom(:) = 0
      count_total      = 0
      num_physdom      = 0
      id_physdom(:)    = 0
      proc_count(:)    = 0
      ncell_offset(:)  = 0
      proc_offset(:)   = 0
      corr_ratio(:)    = 1._wp

      DO j = my_cell_start, my_cell_end
        IF (subset_flag(j) > 0) THEN
          count_physdom(subset_flag(j)) = count_physdom(subset_flag(j)) + 1
          count_total = count_total + 1
        ENDIF
      ENDDO
      IF (p_n_work > 1) THEN
#ifndef NOMPI
        count_physdom(0) = count_total
        CALL mpi_allreduce(mpi_in_place, count_physdom, SIZE(count_physdom), &
             p_int, mpi_sum, p_comm_work, ierror)
        IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
        count_total = count_physdom(0)
#endif
      END IF

      IF (MAXVAL(count_physdom(1:)) < count_total) THEN
        lsplit_merged_domains = .TRUE.
        DO j = 1, max_phys_dom
          IF (count_physdom(j) > 0) THEN
            num_physdom = num_physdom + 1
            id_physdom(num_physdom) = j
            proc_count(num_physdom) = NINT(REAL(count_physdom(j),wp)/REAL(count_total,wp)*REAL(n_proc,wp))
          ENDIF
        ENDDO

        ! Ensure that the sum of partial processor counts matches n_proc
        checksum = SUM(proc_count(1:num_physdom)) - n_proc
        IF (checksum /= 0) THEN
          DO j = 1, num_physdom
            corr_ratio(j) = REAL(proc_count(j),wp) / &
              (REAL(count_physdom(id_physdom(j)),wp)/REAL(count_total,wp)*REAL(n_proc,wp))
          ENDDO
        ENDIF
        IF (checksum > 0) THEN
          DO WHILE (checksum > 0)
            jm = MAXLOC(corr_ratio)
            j = jm(1)
            corr_ratio(j) = 1._wp
            proc_count(j) = proc_count(j) - 1
            checksum = checksum - 1
          ENDDO
        ELSE IF (checksum < 0) THEN
          DO WHILE (checksum < 0)
            jm = MINLOC(corr_ratio)
            j = jm(1)
            corr_ratio(j) = 1._wp
            proc_count(j) = proc_count(j) + 1
            checksum = checksum + 1
          ENDDO
        ENDIF

        ! Compute offset for processor IDs
        DO j = 2, num_physdom
          proc_offset(j) = proc_offset(j-1) + proc_count(j-1)
        ENDDO

        IF(p_pe_work==0) THEN
          WRITE(0,*) 'divide_patch: partial processor counts used for decomposition of merged domain:'
          WRITE(0,*) proc_count(1:num_physdom), 'total: ', n_proc
        ENDIF

      ENDIF

    ENDIF

    IF (.NOT. lsplit_merged_domains) THEN
      num_physdom       = 1
      proc_count(1)     = n_proc
      proc_offset(1)    = 0
      IF (SIZE(subset_flag) > 0) THEN
        id_physdom(1) = p_max(MAXVAL(subset_flag), p_comm_work)
      ELSE
        id_physdom(1) = p_max(0, p_comm_work)
      END IF
    ENDIF


    ! Fill the cell_desc array, it must contain:
    ! cell_desc(1,:)   lat
    ! cell_desc(2,:)   lon
    ! cell_desc(3,:)   cell number (for back-sorting at the end)
    ! cell_desc(4,:)   will be set with the owner
    ALLOCATE(cell_desc(my_cell_start:my_cell_end))

    CALL build_cell_info(my_cell_start, my_cell_end, cell_desc, &
         wrk_p_patch_pre, divide_for_radiation, num_physdom, &
         subset_flag, lsplit_merged_domains, lparent_level, locean, &
         id_physdom, ncell_offset, n_onb_points)
    IF (p_n_work > 1) THEN
#ifndef NOMPI
      CALL mpi_allreduce(ncell_offset, ncell_offset_g, SIZE(ncell_offset), &
           p_int, mpi_sum, wrk_p_patch_pre%dist_array_comm, ierror)
      IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
#endif
    END IF


    DO j = 1, num_physdom

      ncs = ncell_offset(j-1) + my_cell_start
      nce = ncell_offset(j) + my_cell_start - 1
      nc  = ncell_offset(j) - ncell_offset(j-1)

      IF (p_n_work > 1) THEN
#ifndef NOMPI
        nc_g = ncell_offset_g(j) - ncell_offset_g(j-1)
        CALL divide_cells_by_location_mpi(nc_g, nc, cell_desc(ncs:nce), &
             proc_count(j), wrk_p_patch_pre%dist_array_comm)
#endif
      ELSE
        CALL divide_cells_by_location(nc, cell_desc(ncs:nce), &
             0, proc_count(j)-1)
        ! After divide_cells_by_location the cells are sorted by owner,
        ! order them by original cell numbers again
        CALL sort_cell_info_by_cell_number(cell_desc(ncs:nce), nc)
      END IF

      ! Apply shift of processor IDs
      IF (j > 1) THEN
        cell_desc(ncs:nce)%owner = cell_desc(ncs:nce)%owner + proc_offset(j)
      END IF

    ENDDO

    ! set up distributed owner array
    CALL create_dist_cell_owner(dist_cell_owner, &
      &                         wrk_p_patch_pre%n_patch_cells_g, &
      &                         wrk_p_patch_pre%dist_array_comm, &
      &                         wrk_p_patch_pre%dist_array_pes)

    ! Initialization of cell owner field
    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, owner_chunk_ptr)
    owner_chunk_ptr = -1

    IF (p_n_work > 1) THEN
#ifndef NOMPI
      CALL set_owners_mpi(my_cell_start, my_cell_end, dist_cell_owner, &
           wrk_p_patch_pre, cell_desc, ncell_offset, num_physdom, n_onb_points)
#endif
    ELSE
      CALL set_owners(1, wrk_p_patch_pre%n_patch_cells_g, &
           owner_chunk_ptr, wrk_p_patch_pre, cell_desc, ncell_offset, &
           num_physdom, n_onb_points)
    END IF

    DEALLOCATE(cell_desc)
  END SUBROUTINE divide_subset_geometric

  SUBROUTINE build_cell_info(range_start, range_end, cell_desc, &
       wrk_p_patch_pre, divide_for_radiation, num_physdom, subset_flag, &
       lsplit_merged_domains, lparent_level, locean, &
       id_physdom, ncell_offset, n_onb_points)
    INTEGER, INTENT(in) :: range_start, range_end
    TYPE(t_cell_info), INTENT(out) :: cell_desc(range_start:range_end)
    TYPE(t_pre_patch), INTENT(inout) :: wrk_p_patch_pre
    LOGICAL, INTENT(in) :: divide_for_radiation, lsplit_merged_domains, &
         lparent_level, locean
    INTEGER, INTENT(in) :: num_physdom
    ! also ends on range_end, but we want implicit shape argument
    INTEGER, INTENT(in) :: subset_flag(range_start:range_end)
    INTEGER, INTENT(inout) :: ncell_offset(0:max_phys_dom)
    INTEGER, INTENT(out) :: n_onb_points
    INTEGER, INTENT(in) :: id_physdom(max_phys_dom)

    INTEGER :: i, j, nc, jd, idp
#ifdef __BLOCK_GET__
    INTEGER :: k
#endif
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    INTEGER, ALLOCATABLE :: refin_ctrl(:), cell_vertex(:,:)
    REAL(wp), ALLOCATABLE :: temp_lat(:,:), temp_lon(:,:)
    REAL(wp), ALLOCATABLE :: cclat(:), cclon(:)
#else
    INTEGER :: refin_ctrl, cell_vertex
    REAL(wp) :: temp_lat, temp_lon
    REAL(wp) :: cclat, cclon
#endif

    nc = 0
    n_onb_points = 0

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    ALLOCATE(cclat(range_start:range_end), cclon(range_start:range_end), &
      &      cell_vertex(range_start:range_end,3), &
      &      temp_lat(range_start:range_end,3), &
      &      temp_lon(range_start:range_end,3))
#endif

    IF(divide_for_radiation) THEN

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
#ifdef __BLOCK_GET__
      CALL dist_mult_init_blk_comm_dp(2, range_end-range_start+1, 1, wrk_p_patch_pre%cells%dist)
#endif
      DO j = range_start, range_end
        IF (subset_flag(j)<=0) CYCLE ! Cell not in subset
#ifndef __BLOCK_GET__
        CALL dist_mult_array_get( &
          wrk_p_patch_pre%cells%dist, c_center, (/j, 1/), cclat(j))
        CALL dist_mult_array_get(&
          wrk_p_patch_pre%cells%dist, c_center, (/j, 2/), cclon(j))
#else
        i = j - range_start + 1
        CALL dist_mult_array_get_dp_blk( &
             wrk_p_patch_pre%cells%dist, c_center, (/j, 1/), cclat(j), &
             proc_lst_dp(1,i,1), ind_lst_dp(1,i,1))
        CALL dist_mult_array_get_dp_blk( &
             wrk_p_patch_pre%cells%dist, c_center, (/j, 2/), cclon(j), &
             proc_lst_dp(2,i,1), ind_lst_dp(2,i,1))
#endif
      END DO
#ifdef __BLOCK_GET__
      call dist_mult_do_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)
#ifdef __BLOCK_GET__
      DO j = range_start, range_end
        IF (subset_flag(j)<=0) CYCLE ! Cell not in subset
        i = j - range_start + 1
        IF (proc_lst_dp(1,i,1) /= -1) THEN
          cclat(j) = recv_array_dp(proc_lst_dp(1,i,1))%lst(ind_lst_dp(1,i,1))
        END IF
        IF (proc_lst_dp(2,i,1) /= -1) THEN
          cclon(j) = recv_array_dp(proc_lst_dp(2,i,1))%lst(ind_lst_dp(2,i,1))
        END IF
      END DO
      CALL dist_mult_end_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif
      DO j = range_start, range_end
        cell_desc(j)%lat = HUGE(cell_desc(j)%lat)! for finding min lat/lon
        cell_desc(j)%lon = HUGE(cell_desc(j)%lon) ! for finding min lat/lon
        IF (subset_flag(j)<=0) CYCLE ! Cell not in subset

        ! Patch division for radiation calculations:
        ! To minimize load imbalance, every patch contains 10 areas
        ! distributed in a way similar as the "diamonds" in GME
        ! This is accomplished by mapping all cells to one section
        ! lying in the NH and having a width of 0.4*pi (72 deg)

        IF (cclat(j)>=0._wp .AND. cclon(j)>=-0.2_wp*pi .AND. cclon(j)<=0.2_wp*pi) THEN
        ELSE IF (cclat(j)>=0._wp .AND. cclon(j)>=0.2_wp*pi .AND. cclon(j)<=0.6_wp*pi) THEN
          cclon(j) = cclon(j) - 0.4_wp*pi
        ELSE IF (cclat(j)>=0._wp .AND. cclon(j)>=0.6_wp*pi) THEN
          cclon(j) = cclon(j) - 0.8_wp*pi
        ELSE IF (cclat(j)>=0._wp .AND. cclon(j)<=-0.6_wp*pi) THEN
          cclon(j) = cclon(j) + 0.8_wp*pi
        ELSE IF (cclat(j)>=0._wp .AND. cclon(j)>=-0.6_wp*pi .AND. &
          &      cclon(j)<=-0.2_wp*pi) THEN
          cclon(j) = cclon(j) + 0.4_wp*pi
        ELSE IF (cclat(j)<0._wp .AND. &
          &      (cclon(j)<=-0.8_wp*pi .OR. cclon(j)>=0.8_wp*pi)) THEN
          cclat(j) = -cclat(j)
          cclon(j) = cclon(j) + pi
        ELSE IF (cclat(j)<0._wp .AND. cclon(j)>=-0.8_wp*pi .AND. &
          &      cclon(j)<=-0.4_wp*pi) THEN
          cclat(j) = -cclat(j)
          cclon(j) = cclon(j) + 0.6_wp*pi
        ELSE IF (cclat(j)<0._wp .AND. cclon(j)>=-0.4_wp*pi .AND. &
          &      cclon(j)<=0.0_wp*pi) THEN
          cclat(j) = -cclat(j)
          cclon(j) = cclon(j) + 0.2_wp*pi
        ELSE IF (cclat(j)<0._wp .AND. cclon(j)>=0.0_wp*pi .AND. &
          &      cclon(j)<=0.4_wp*pi) THEN
          cclat(j) = -cclat(j)
          cclon(j) = cclon(j) - 0.2_wp*pi
        ELSE IF (cclat(j)<0._wp .AND. cclon(j)>=0.4_wp*pi .AND. &
          &      cclon(j)<=0.8_wp*pi) THEN
          cclat(j) = -cclat(j)
          cclon(j) = cclon(j) - 0.6_wp*pi
        ENDIF

        IF (cclon(j) > pi) THEN
          cclon(j) = cclon(j) - 2._wp*pi
          cclon(j) = cclon(j) + 2._wp*pi
        ELSE IF (cclon(j) < -pi) THEN
        ENDIF

        cell_desc(range_start + nc)%lat = fxp_lat(cclat(j))
        cell_desc(range_start + nc)%lon = fxp_lon(cclon(j))
        cell_desc(range_start + nc)%cell_number = j
        cell_desc(range_start + nc)%owner = 0

        nc = nc+1 ! Cell counter

      ENDDO
#else
      DO j = range_start, range_end
        cell_desc(j)%lat = HUGE(cell_desc(j)%lat)! for finding min lat/lon
        cell_desc(j)%lon = HUGE(cell_desc(j)%lon) ! for finding min lat/lon
        IF (subset_flag(j)<=0) CYCLE ! Cell not in subset

        ! Patch division for radiation calculations:
        ! To minimize load imbalance, every patch contains 10 areas
        ! distributed in a way similar as the "diamonds" in GME
        ! This is accomplished by mapping all cells to one section
        ! lying in the NH and having a width of 0.4*pi (72 deg)

        CALL dist_mult_array_get( &
          wrk_p_patch_pre%cells%dist, c_center, (/j, 1/), cclat)
        CALL dist_mult_array_get(&
          wrk_p_patch_pre%cells%dist, c_center, (/j, 2/), cclon)

        IF (cclat>=0._wp .AND. cclon>=-0.2_wp*pi .AND. cclon<=0.2_wp*pi) THEN
        ELSE IF (cclat>=0._wp .AND. cclon>=0.2_wp*pi .AND. cclon<=0.6_wp*pi) THEN
          cclon = cclon - 0.4_wp*pi
        ELSE IF (cclat>=0._wp .AND. cclon>=0.6_wp*pi) THEN
          cclon = cclon - 0.8_wp*pi
        ELSE IF (cclat>=0._wp .AND. cclon<=-0.6_wp*pi) THEN
          cclon = cclon + 0.8_wp*pi
        ELSE IF (cclat>=0._wp .AND. cclon>=-0.6_wp*pi .AND. cclon<=-0.2_wp*pi) THEN
          cclon = cclon + 0.4_wp*pi
        ELSE IF (cclat<0._wp .AND. (cclon<=-0.8_wp*pi .OR. cclon>=0.8_wp*pi)) THEN
          cclat = -cclat
          cclon = cclon + pi
        ELSE IF (cclat<0._wp .AND. cclon>=-0.8_wp*pi .AND. cclon<=-0.4_wp*pi) THEN
          cclat = -cclat
          cclon = cclon + 0.6_wp*pi
        ELSE IF (cclat<0._wp .AND. cclon>=-0.4_wp*pi .AND. cclon<=0.0_wp*pi) THEN
          cclat = -cclat
          cclon = cclon + 0.2_wp*pi
        ELSE IF (cclat<0._wp .AND. cclon>=0.0_wp*pi .AND. cclon<=0.4_wp*pi) THEN
          cclat = -cclat
          cclon = cclon - 0.2_wp*pi
        ELSE IF (cclat<0._wp .AND. cclon>=0.4_wp*pi .AND. cclon<=0.8_wp*pi) THEN
          cclat = -cclat
          cclon = cclon - 0.6_wp*pi
        ENDIF

        IF (cclon > pi) THEN
          cclon = cclon - 2._wp*pi
          cclon = cclon + 2._wp*pi
        ELSE IF (cclon < -pi) THEN
        ENDIF

        cell_desc(range_start + nc)%lat = fxp_lat(cclat)
        cell_desc(range_start + nc)%lon = fxp_lon(cclon)
        cell_desc(range_start + nc)%cell_number = j
        cell_desc(range_start + nc)%owner = 0

        nc = nc+1 ! Cell counter

      ENDDO
#endif

      ncell_offset(0) = 0
      ncell_offset(1) = nc

    ELSE IF (wrk_p_patch_pre%cells%max_connectivity == 6) THEN

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
#ifdef __BLOCK_GET__
      CALL dist_mult_init_blk_comm_dp(2, range_end-range_start+1, 1, wrk_p_patch_pre%cells%dist)
#endif
      DO j = range_start, range_end

        IF (subset_flag(j) <= 0) CYCLE

#ifndef __BLOCK_GET__
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
             (/j, 1/), cclat(j))
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
             (/j, 2/), cclon(j))
#else
        i = j - range_start + 1
        CALL dist_mult_array_get_dp_blk( &
             wrk_p_patch_pre%cells%dist, c_center, (/j, 1/), cclat(j), &
             proc_lst_dp(1,i,1), ind_lst_dp(1,i,1))
        CALL dist_mult_array_get_dp_blk( &
             wrk_p_patch_pre%cells%dist, c_center, (/j, 2/), cclon(j), &
             proc_lst_dp(2,i,1), ind_lst_dp(2,i,1))
#endif
      END DO
#ifdef __BLOCK_GET__
      call dist_mult_do_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif
      CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)
#ifdef __BLOCK_GET__
      DO j = range_start, range_end
        IF (subset_flag(j)<=0) CYCLE ! Cell not in subset
        i = j - range_start + 1
        IF (proc_lst_dp(1,i,1) /= -1) THEN
          cclat(j) = recv_array_dp(proc_lst_dp(1,i,1))%lst(ind_lst_dp(1,i,1))
        END IF
        IF (proc_lst_dp(2,i,1) /= -1) THEN
          cclon(j) = recv_array_dp(proc_lst_dp(2,i,1))%lst(ind_lst_dp(2,i,1))
        END IF
      END DO
      CALL dist_mult_end_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif
      DO j = range_start, range_end

        IF (subset_flag(j) <= 0) CYCLE
        cell_desc(range_start + nc)%lat = fxp_lat(cclat(j))
        cell_desc(range_start + nc)%lon = fxp_lon(cclon(j))
        cell_desc(range_start + nc)%cell_number = j
        cell_desc(range_start + nc)%owner = 0

        nc = nc + 1 ! Cell counter
      ENDDO
#else
      DO j = range_start, range_end

        IF (subset_flag(j) <= 0) CYCLE

        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
             (/j, 1/), cclat)
        CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
             (/j, 2/), cclon)
        cell_desc(range_start + nc)%lat = fxp_lat(cclat)
        cell_desc(range_start + nc)%lon = fxp_lon(cclon)
        cell_desc(range_start + nc)%cell_number = j
        cell_desc(range_start + nc)%owner = 0

        nc = nc + 1 ! Cell counter
      ENDDO
#endif

      ncell_offset(1) = nc

    ELSE ! domain decomposition for triangular cells with optional splitting into physical domains

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      ALLOCATE(refin_ctrl(range_start:range_end))
#endif

      DO jd = 1, num_physdom

        idp = id_physdom(jd)

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
#ifdef __BLOCK_GET__
        CALL dist_mult_init_blk_comm(1, range_end-range_start+1, 1, wrk_p_patch_pre%cells%dist)
#endif
        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_refin_ctrl, &
               (/j/), refin_ctrl(j))
#else
          i = j - range_start + 1
          CALL dist_mult_array_get_i4_blk( &
               wrk_p_patch_pre%cells%dist, c_refin_ctrl, (/j/), refin_ctrl(j), &
               proc_lst(1,i,1), ind_lst(1,i,1))
#endif

        END DO
#ifdef __BLOCK_GET__
        CALL dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
#endif
        CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)
#ifdef __BLOCK_GET__
      DO j = range_start, range_end
        IF (subset_flag(j)<=0) CYCLE ! Cell not in subset
        i = j - range_start + 1
        IF (proc_lst(1,i,1) /= -1) THEN
          refin_ctrl(j) = recv_array(proc_lst(1,i,1))%lst(ind_lst(1,i,1))
        END IF
      END DO
      CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
#endif

#ifdef __BLOCK_GET__
        call dist_mult_init_blk_comm_dp(2, range_end-range_start+1, 1, &
             wrk_p_patch_pre%cells%dist)
        call dist_mult_init_blk_comm(3, range_end-range_start+1, 1, &
             wrk_p_patch_pre%cells%dist)
#endif
        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

          IF (lparent_level .AND. refin_ctrl(j) == -1 .OR. .NOT. lparent_level &
            & .AND. refin_ctrl(j) >= 1 .AND. refin_ctrl(j) <= 3 .AND. .NOT. locean) &
            CYCLE

#ifndef __BLOCK_GET__
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
               (/j, 1/), cclat(j))
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
               (/j, 2/), cclon(j))

          DO i = 1, 3
            CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                 (/j,i/), cell_vertex(j,i))
          ENDDO
#else
          k = j - range_start + 1
          CALL dist_mult_array_get_dp_blk( &
             wrk_p_patch_pre%cells%dist, c_center, (/j, 1/), cclat(j), &
             proc_lst_dp(1,k,1), ind_lst_dp(1,k,1))
          CALL dist_mult_array_get_dp_blk( &
             wrk_p_patch_pre%cells%dist, c_center, (/j, 2/), cclon(j), &
             proc_lst_dp(2,k,1), ind_lst_dp(2,k,1))
          DO i = 1, 3
            CALL dist_mult_array_get_i4_blk( &
               wrk_p_patch_pre%cells%dist, c_vertex, (/j,i /), cell_vertex(j,i), &
               proc_lst(i,k,1), ind_lst(i,k,1))
          END DO
#endif
        END DO
#ifdef __BLOCK_GET__
        call dist_mult_do_blk_comm(wrk_p_patch_pre%cells%dist)
        call dist_mult_do_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif
        CALL dist_mult_array_rma_sync(wrk_p_patch_pre%cells%dist)

#ifdef __BLOCK_GET__
        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current
          ! physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE
          IF (lparent_level .AND. refin_ctrl(j) == -1 .OR. .NOT. lparent_level &
            & .AND. refin_ctrl(j) >= 1 .AND. refin_ctrl(j) <= 3 .AND. .NOT. locean) CYCLE
          k = j - range_start + 1
          !doubles
          IF (proc_lst_dp(1,k,1) /= -1) THEN
            cclat(j) = recv_array_dp(proc_lst_dp(1,k,1))%lst(ind_lst_dp(1,k,1))
          END IF
          IF (proc_lst_dp(2,k,1) /= -1) THEN
            cclon(j) = recv_array_dp(proc_lst_dp(2,k,1))%lst(ind_lst_dp(2,k,1))
          END IF
          !ints
          DO i = 1, 3
            IF (proc_lst(i,k,1) /= -1) THEN
              cell_vertex(j,i) = recv_array(proc_lst(i,k,1))%lst(ind_lst(i,k,1))
            END IF
          END DO
        END DO
        CALL dist_mult_end_blk_comm(wrk_p_patch_pre%cells%dist)
        CALL dist_mult_end_blk_comm_dp(wrk_p_patch_pre%cells%dist)
#endif
#ifdef __BLOCK_GET__
        call dist_mult_init_blk_comm_dp(2, 3, range_end-range_start+1, &
             wrk_p_patch_pre%verts%dist)
#endif
        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

          IF (lparent_level .AND. refin_ctrl(j) == -1 .OR. .NOT. lparent_level &
            & .AND. refin_ctrl(j) >= 1 .AND. refin_ctrl(j) <= 3 .AND. .NOT. locean) &
            CYCLE
#ifdef __BLOCK_GET__
          k = j - range_start + 1
#endif

          DO i = 1, 3
#ifndef __BLOCK_GET__
            CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
                 (/cell_vertex(j,i), 1/), temp_lat(j,i))
            CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
                 (/cell_vertex(j,i), 2/), temp_lon(j,i))
#else
            CALL dist_mult_array_get_dp_blk(wrk_p_patch_pre%verts%dist, &
                 v_vertex, (/cell_vertex(j,i), 1/), temp_lat(j,i), &
                 proc_lst_dp(1,i,k), ind_lst_dp(1,i,k))
            CALL dist_mult_array_get_dp_blk(wrk_p_patch_pre%verts%dist, &
                 v_vertex, (/cell_vertex(j,i), 2/), temp_lon(j,i), &
                 proc_lst_dp(2,i,k), ind_lst_dp(2,i,k))
#endif
          ENDDO
        END DO
#ifdef __BLOCK_GET__
        call dist_mult_do_blk_comm_dp(wrk_p_patch_pre%verts%dist)
#endif
        CALL dist_mult_array_rma_sync(wrk_p_patch_pre%verts%dist)
#ifdef __BLOCK_GET__

        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current
          ! physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

          IF (lparent_level .AND. refin_ctrl(j) == -1 .OR. .NOT. lparent_level &
            & .AND. refin_ctrl(j) >= 1 .AND. refin_ctrl(j) <= 3 .AND. .NOT.locean) &
            CYCLE
          k = j - range_start + 1
          DO i = 1, 3
            IF (proc_lst_dp(1,i,k) /= -1) THEN
              temp_lat(j,i) = recv_array_dp(proc_lst_dp(1,i,k))%lst(ind_lst_dp(1,i,k))
            END IF
            IF (proc_lst_dp(2,i,k) /= -1) THEN
              temp_lon(j,i) = recv_array_dp(proc_lst_dp(2,i,k))%lst(ind_lst_dp(2,i,k))
            END IF
          END DO
        END DO
        call dist_mult_end_blk_comm_dp(wrk_p_patch_pre%verts%dist)
#endif

!CDIR NODEP
        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

          ! Disregard outer nest boundary points for the time being. They do very little
          ! computational work, so they can be added to the closest PEs afterwards
          IF (lparent_level .AND. refin_ctrl(j) == -1 .OR. .NOT. lparent_level &
            & .AND. refin_ctrl(j) >= 1 .AND. refin_ctrl(j) <= 3 .AND. .NOT. locean) THEN
            cell_desc(range_end-n_onb_points)%cell_number = j
            n_onb_points = n_onb_points + 1
            CYCLE
          ENDIF

          ! Using the center of the cells for geometric subdivision leads
          ! to "toothed" edges of the subdivision area
          ! Thus we use the minimum lat/lon as subdivision criterion.

          IF (cclat(j) >= 0._wp) THEN
            DO i = 1, 3
              cclat(j) = MAX(cclat(j), temp_lat(j,i))
              cclon(j) = MAX(cclon(j), temp_lon(j,i))
            ENDDO
          ELSE
            DO i = 1, 3
              cclat(j) = MIN(cclat(j), temp_lat(j,i))
              cclon(j) = MAX(cclon(j), temp_lon(j,i))
            ENDDO
          ENDIF

          cell_desc(range_start + nc)%lat = fxp_lat(cclat(j))
          cell_desc(range_start + nc)%lon = fxp_lon(cclon(j))
          cell_desc(range_start + nc)%cell_number = j
          cell_desc(range_start + nc)%owner = 0

          nc = nc+1 ! Cell counter

        ENDDO
#else

!CDIR NODEP
        DO j = range_start, range_end

          ! Skip cell if it is not in subset or does not belong to current physical domain
          IF (subset_flag(j) /= idp .AND. lsplit_merged_domains .OR. subset_flag(j) <= 0) CYCLE

          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_refin_ctrl, &
               (/j/), refin_ctrl)

          ! Disregard outer nest boundary points for the time being. They do very little
          ! computational work, so they can be added to the closest PEs afterwards
          IF (lparent_level .AND. refin_ctrl == -1 .OR. .NOT. lparent_level &
            & .AND. refin_ctrl >= 1 .AND. refin_ctrl <= 3 .AND. .NOT. locean) THEN
            cell_desc(range_end-n_onb_points)%cell_number = j
            n_onb_points = n_onb_points + 1
            CYCLE
          ENDIF

          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
               (/j, 1/), cclat)
          CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_center, &
               (/j, 2/), cclon)

          ! Using the center of the cells for geometric subdivision leads
          ! to "toothed" edges of the subdivision area
          ! Thus we use the minimum lat/lon as subdivision criterion.

          IF (cclat >= 0._wp) THEN
            DO i = 1, 3
              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/j,i/), cell_vertex)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
                   (/cell_vertex, 1/), temp_lat)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
                   (/cell_vertex, 2/), temp_lon)
              cclat = MAX(cclat, temp_lat)
              cclon = MAX(cclon, temp_lon)
            ENDDO
          ELSE
            DO i = 1, 3
              CALL dist_mult_array_get(wrk_p_patch_pre%cells%dist, c_vertex, &
                   (/j,i/), cell_vertex)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
                   (/cell_vertex, 1/), temp_lat)
              CALL dist_mult_array_get(wrk_p_patch_pre%verts%dist, v_vertex, &
                   (/cell_vertex, 2/), temp_lon)
              cclat = MIN(cclat, temp_lat)
              cclon = MAX(cclon, temp_lon)
            ENDDO
          ENDIF

          cell_desc(range_start + nc)%lat = fxp_lat(cclat)
          cell_desc(range_start + nc)%lon = fxp_lon(cclon)
          cell_desc(range_start + nc)%cell_number = j
          cell_desc(range_start + nc)%owner = 0

          nc = nc+1 ! Cell counter

        ENDDO

#endif

        ncell_offset(jd) = nc

      ENDDO

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      DEALLOCATE(refin_ctrl)
#endif

    ENDIF

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    DEALLOCATE(cclat, cclon, cell_vertex, temp_lat, temp_lon)
#endif

  END SUBROUTINE build_cell_info

  ! Set owner list (of complete patch)
  SUBROUTINE set_owners(range_start, range_end, &
       owner, wrk_p_patch_pre, cell_desc, ncell_offset, &
       num_physdom, n_onb_points)
    INTEGER, INTENT(in) :: range_start, range_end
     !> receives the owner PE for every cell
    INTEGER, INTENT(out) :: owner(range_start:range_end)
    TYPE(t_pre_patch), INTENT(inout) :: wrk_p_patch_pre
    TYPE(t_cell_info), INTENT(inout) :: cell_desc(range_start:range_end)
    !> if > 0 a cell belongs to the subset
    INTEGER, INTENT(in) :: ncell_offset(0:max_phys_dom)
    INTEGER, INTENT(in) :: num_physdom, n_onb_points

    INTEGER :: i, ii, j, jd, jj, jn, &
         first_unassigned_onb, last_unassigned_onb, &
         ncs, nce, num_set_on_this_pass
    TYPE(t_cell_info) :: temp
    INTEGER, POINTER :: cells_num_edges_local_ptr(:), &
      &                 cells_neighbor_local_ptr(:,:)
    INTEGER, ALLOCATABLE :: set_on_this_pass(:)

    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_num_edges, &
      &                            cells_num_edges_local_ptr)
    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_neighbor, &
      &                            cells_neighbor_local_ptr)

    DO jd = 1, num_physdom
      ncs = ncell_offset(jd-1)+1
      nce = ncell_offset(jd)
      DO j = ncs, nce
        owner(cell_desc(j)%cell_number) = cell_desc(j)%owner
      END DO
    ENDDO

    ! Add outer nest boundary points that have been disregarded so far
    IF (n_onb_points > 0) THEN
      first_unassigned_onb = range_end - n_onb_points + 1
      last_unassigned_onb = range_end
      ! Iterations are needed because outer nest boundary row contains
      ! indirect neighbors
      ALLOCATE(set_on_this_pass(n_onb_points))
      DO WHILE (first_unassigned_onb <= last_unassigned_onb)
        num_set_on_this_pass = 0
        DO i = last_unassigned_onb, first_unassigned_onb, -1
          j = cell_desc(i)%cell_number
          DO ii = 1, cells_num_edges_local_ptr(j)
            jn = cells_neighbor_local_ptr(j, ii)
            IF (jn > 0) THEN
              IF (owner(jn) >= 0) THEN
                ! move found cells to end of onb area
                temp = cell_desc(i)
                DO jj = i, last_unassigned_onb - 1
                  cell_desc(jj) = cell_desc(jj + 1)
                END DO
                cell_desc(last_unassigned_onb) = temp
                owner(j) = -owner(jn) - 1
                num_set_on_this_pass = num_set_on_this_pass + 1
                set_on_this_pass(num_set_on_this_pass) = j
                last_unassigned_onb = last_unassigned_onb - 1
                EXIT
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO i = 1, num_set_on_this_pass
          j = set_on_this_pass(i)
          owner(j) = -owner(j) - 1
        END DO
      ENDDO
    ENDIF

  END SUBROUTINE set_owners

#ifndef NOMPI
  ! Set owner list (of complete patch)
  SUBROUTINE set_owners_mpi(range_start, range_end, dist_cell_owner, &
       wrk_p_patch_pre, cell_desc, ncell_offset, num_physdom, n_onb_points)
    USE mpi
    USE ppm_extents
    USE ppm_distributed_array
    INTEGER, INTENT(in) :: range_start, range_end
     !> receives the owner PE for every cell
    TYPE(dist_mult_array), INTENT(inout) :: dist_cell_owner
    TYPE(t_pre_patch), INTENT(inout) :: wrk_p_patch_pre
    TYPE(t_cell_info), INTENT(inout) :: cell_desc(range_start:range_end)
    !> if > 0 a cell belongs to the subset
    INTEGER, INTENT(in) :: ncell_offset(0:max_phys_dom)
    INTEGER, INTENT(in) :: num_physdom, n_onb_points

    INTEGER :: i, ii, j, jd, jj, jn, first_unassigned_onb, &
      &        last_unassigned_onb, ncs, nce, owner_idx(1), ierror
    TYPE(t_cell_info) :: temp
    INTEGER, POINTER :: owners_local_ptr(:), cells_num_edges_local_ptr(:), &
      &                 cells_neighbor_local_ptr(:,:)
    INTEGER, ALLOCATABLE :: tmp_owners_local(:)
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    INTEGER, ALLOCATABLE :: owner_value(:,:)
#else
    INTEGER :: owner_value
#endif
    LOGICAL :: owner_missing, any_onb_points
    CHARACTER(14), PARAMETER :: routine = 'set_owners_mpi'

    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, owners_local_ptr)
    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_num_edges, &
      &                            cells_num_edges_local_ptr)
    CALL dist_mult_array_local_ptr(wrk_p_patch_pre%cells%dist, c_neighbor, &
      &                            cells_neighbor_local_ptr)

    DO jd = 1, num_physdom
      ncs = ncell_offset(jd-1) + range_start
      nce = ncell_offset(jd) + range_start - 1
      DO j = ncs, nce
        owners_local_ptr(cell_desc(j)%cell_number) = cell_desc(j)%owner
      END DO
    END DO

    ALLOCATE(tmp_owners_local(range_start:range_end))
    tmp_owners_local = owners_local_ptr

    ! Add outer nest boundary points that have been disregarded so far
    any_onb_points = n_onb_points > 0
    CALL mpi_allreduce(mpi_in_place, any_onb_points, 1, p_bool, &
         mpi_lor, wrk_p_patch_pre%dist_array_comm, ierror)
    IF (any_onb_points) THEN

      first_unassigned_onb = range_end - n_onb_points + 1
      last_unassigned_onb = range_end

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      ALLOCATE(owner_value(first_unassigned_onb:last_unassigned_onb, &
        &                  wrk_p_patch_pre%cells%max_connectivity))
      ! Iterations are needed because outer nest boundary row contains
      ! indirect neighbors
      owner_missing = first_unassigned_onb <= last_unassigned_onb
      CALL mpi_allreduce(mpi_in_place, owner_missing, 1, p_bool, &
           mpi_lor, wrk_p_patch_pre%dist_array_comm, ierror)
      DO WHILE (owner_missing)
        IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
        CALL dist_mult_array_expose(dist_cell_owner)
        DO i = last_unassigned_onb, first_unassigned_onb, -1
          j = cell_desc(i)%cell_number
          DO ii = 1, cells_num_edges_local_ptr(j)
            jn = cells_neighbor_local_ptr(j, ii)
            IF (jn > 0) THEN
              owner_idx(1) = jn
              CALL dist_mult_array_get(dist_cell_owner, 1, owner_idx, &
                &                      owner_value(i,ii))
            ENDIF
          ENDDO
        ENDDO
        CALL dist_mult_array_rma_sync(dist_cell_owner)
        DO i = last_unassigned_onb, first_unassigned_onb, -1
          j = cell_desc(i)%cell_number
          DO ii = 1, cells_num_edges_local_ptr(j)
            jn = cells_neighbor_local_ptr(j, ii)
            IF (jn > 0) THEN
              IF (owner_value(i,ii) >= 0) THEN
                ! move found cells to end of onb area
                temp = cell_desc(i)
                DO jj = i, last_unassigned_onb - 1
                  cell_desc(jj) = cell_desc(jj + 1)
                END DO
                cell_desc(last_unassigned_onb) = temp
                tmp_owners_local(j) = owner_value(i,ii)
                last_unassigned_onb = last_unassigned_onb - 1
                EXIT
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        CALL dist_mult_array_unexpose(dist_cell_owner)
        owners_local_ptr = tmp_owners_local
        owner_missing = first_unassigned_onb <= last_unassigned_onb
        CALL mpi_allreduce(mpi_in_place, owner_missing, 1, p_bool, &
             mpi_lor, wrk_p_patch_pre%dist_array_comm, ierror)
      ENDDO
      DEALLOCATE(owner_value)
#else

      ! Iterations are needed because outer nest boundary row contains
      ! indirect neighbours
      owner_missing = first_unassigned_onb <= last_unassigned_onb
      CALL mpi_allreduce(mpi_in_place, owner_missing, 1, p_bool, &
           mpi_lor, wrk_p_patch_pre%dist_array_comm, ierror)
      DO WHILE (owner_missing)
        IF (ierror /= mpi_success) CALL finish(routine, 'error in allreduce')
        CALL dist_mult_array_expose(dist_cell_owner)
        DO i = last_unassigned_onb, first_unassigned_onb, -1
          j = cell_desc(i)%cell_number
          DO ii = 1, cells_num_edges_local_ptr(j)
            jn = cells_neighbor_local_ptr(j, ii)
            IF (jn > 0) THEN
              owner_idx(1) = jn
              CALL dist_mult_array_get(dist_cell_owner, 1, owner_idx, &
                &                      owner_value)
              IF (owner_value >= 0) THEN
                ! move found cells to end of onb area
                temp = cell_desc(i)
                DO jj = i, last_unassigned_onb - 1
                  cell_desc(jj) = cell_desc(jj + 1)
                END DO
                cell_desc(last_unassigned_onb) = temp
                tmp_owners_local(j) = owner_value
                last_unassigned_onb = last_unassigned_onb - 1
                EXIT
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        CALL dist_mult_array_unexpose(dist_cell_owner)
        owners_local_ptr = tmp_owners_local
        owner_missing = first_unassigned_onb <= last_unassigned_onb
        CALL mpi_allreduce(mpi_in_place, owner_missing, 1, p_bool, &
             mpi_lor, wrk_p_patch_pre%dist_array_comm, ierror)
      ENDDO

#endif

    ENDIF

  END SUBROUTINE set_owners_mpi
#endif
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  !>
  SUBROUTINE read_netcdf_decomposition(netcdf_file_name, dist_cell_owner, &
    &                                  no_of_cells, dist_array_comm, &
    &                                  dist_array_pes)
    CHARACTER(LEN=*) :: netcdf_file_name
    TYPE(dist_mult_array), INTENT(OUT) :: dist_cell_owner
    INTEGER, INTENT(IN) :: no_of_cells, dist_array_comm
    TYPE(extent), INTENT(IN) :: dist_array_pes

    TYPE(global_array_desc) :: dist_cell_owner_desc(1)
    TYPE(extent) :: local_chunk(1,1)
    INTEGER, POINTER :: local_ptr(:)
    INTEGER, ALLOCATABLE, TARGET :: buffer(:,:)
    INTEGER :: i, ncid
    TYPE(t_grid_domain_decomp_info) :: decomp_info
    TYPE(t_distrib_read_data) :: io_data
    INTEGER :: dist_array_pes_start, dist_array_pes_size
    CHARACTER(*), PARAMETER :: method_name = "read_netcdf_decomposition"
    TYPE(t_ptr_2d_int) :: tmp_ptr(1)

    dist_cell_owner_desc(1)%a_rank = 1
    dist_cell_owner_desc(1)%rect(1)%first = 1
    dist_cell_owner_desc(1)%rect(1)%size = no_of_cells
    dist_cell_owner_desc(1)%element_dt = ppm_int

    dist_array_pes_start = extent_start(dist_array_pes) - 1
    dist_array_pes_size = extent_size(dist_array_pes)

    ! compute uniform partition
    local_chunk(1,1) = uniform_partition(dist_cell_owner_desc(1)%rect(1), &
        &                                dist_array_pes_size, &
        &                                p_pe_work + 1 - dist_array_pes_start)

    ! generate distributed array cell_owner array
    dist_cell_owner = dist_mult_array_new(dist_cell_owner_desc, local_chunk, &
      &          dist_array_comm, &
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
      &          sync_mode=sync_mode_active_target &
#else
      &          cache_size=MIN(10, CEILING(SQRT(REAL(dist_array_pes_size)))) &
#endif
    )

    ! get pointer to local chunk of cell_owner array
    CALL dist_mult_array_local_ptr(dist_cell_owner, 1, local_ptr)

    ! allocate temporary buffer (read method only outputs nproma-blocked data)
    ALLOCATE(buffer(nproma, (SIZE(local_ptr) + nproma - 1)/nproma))

    ! setup io data
    ALLOCATE(decomp_info%glb_index(SIZE(local_ptr)))
    decomp_info%glb_index = (/(i, i=LBOUND(local_ptr,1), UBOUND(local_ptr,1))/)
    CALL setup_distrib_read(no_of_cells, decomp_info, io_data)

    ! read cell owner
    IF (p_pe_work == 0) &
      WRITE(0,*) "Read decomposition from file: ", TRIM(netcdf_file_name)
    ncid = distrib_nf_open(netcdf_file_name)
    tmp_ptr(1)%p => buffer
    CALL distrib_read(ncid, "cell_owner", tmp_ptr, (/io_data/))
    CALL distrib_nf_close(ncid)

    DO i = 0, SIZE(local_ptr) - 1
      local_ptr(i + LBOUND(local_ptr, 1)) = &
        buffer(MOD(i, nproma) + 1, i / nproma + 1)
    END DO

    CALL dist_mult_array_expose(dist_cell_owner)

    ! clean-up
    DEALLOCATE(decomp_info%glb_index)
    CALL delete_distrib_read(io_data)

  END SUBROUTINE read_netcdf_decomposition
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE write_netcdf_decomposition(netcdf_file_name, dist_cell_owner, &
    &                                   no_of_cells)
    CHARACTER(*) :: netcdf_file_name
    TYPE(dist_mult_array), INTENT(INOUT) :: dist_cell_owner
    INTEGER, INTENT(IN) :: no_of_cells

    INTEGER :: i, j, ncid, dimid, varid, part_size, num_parts
    INTEGER :: curr_part_size, start_value(1), value_count(1)
    INTEGER, ALLOCATABLE :: cell_owner_buffer(:)
    CHARACTER(*), PARAMETER :: routine = "write_netcdf_decomposition"

    IF (p_pe_work == 0) THEN
      WRITE(0,*) "Write decomposition to file: ", TRIM(netcdf_file_name)
      ! generate decomposition file
      CALL nf(nf_create(TRIM(netcdf_file_name), NF_CLOBBER, ncid), routine)

      ! define dimension
      CALL nf(nf_def_dim(ncid, "ncells", no_of_cells, dimid), routine)

      ! define variable
      CALL nf(nf_def_var(ncid, "cell_owner", NF_INT, 1, (/dimid/), varid), routine)

      ! put the number of processes
      CALL nf(nf_put_att_int(ncid, varid, "nprocs", NF_INT, 1, (/p_n_work/)), routine)

      CALL nf(nf_enddef(ncid), routine)

      ! write cell owner to file
      part_size = MIN(1024 * 1024 / 8, no_of_cells);
      num_parts = (no_of_cells + part_size - 1) / part_size
      ALLOCATE(cell_owner_buffer(part_size))
      DO i = 1, num_parts
        curr_part_size = MIN(part_size, no_of_cells - part_size * (i - 1))
        DO j = 1, curr_part_size
          CALL dist_mult_array_get(dist_cell_owner, 1, &
            &                      (/j + part_size * (i - 1)/), &
            &                      cell_owner_buffer(j))
        END DO
#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
        CALL dist_mult_array_rma_sync(dist_cell_owner)
#endif
        start_value(1) = part_size * (i - 1) + 1
        value_count(1) = curr_part_size
        CALL nf(nf_put_vara_int(ncid, varid, start_value, value_count, &
          &     cell_owner_buffer), routine)
      END DO
      DEALLOCATE(cell_owner_buffer)

      CALL nf(nf_close(ncid), routine)

#ifdef HAVE_SLOW_PASSIVE_TARGET_ONESIDED
    ELSE
      part_size = MIN(1024 * 1024 / 8, no_of_cells);
      num_parts = (no_of_cells + part_size - 1) / part_size

      DO i = 1, num_parts
        CALL dist_mult_array_rma_sync(dist_cell_owner)
      END DO
#endif

    END IF

  END SUBROUTINE write_netcdf_decomposition
  !-------------------------------------------------------------------------

END MODULE mo_setup_subdivision
