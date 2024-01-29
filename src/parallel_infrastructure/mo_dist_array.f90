#ifndef NOMPI
! Distributed array implementation, optimized for read-only access in
! parallel phase
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

#if defined (__SX__) || defined (__NEC_VH__)
#define __SXorVH__
#endif

MODULE ppm_distributed_array
  USE mo_kind, ONLY: i4, i8, sp, dp
  USE mo_mpi, ONLY: mp_i4_extent => p_int_i4_byte, &
       mp_i8_extent => p_int_i8_byte, &
       mp_l_extent => p_int_byte, &
       mp_sp_extent => p_real_sp_byte, &
       mp_dp_extent => p_real_dp_byte, &
       ppm_real_dp => p_real_dp, &
       ppm_real_sp => p_real_sp, &
       ppm_int => p_int, &
       ppm_int_i8 => p_int_i8, &
       ppm_bool => p_bool, p_pe
#ifdef __BLOCK_GET__
  USE mo_mpi, ONLY: p_comm_work
#endif
  USE ppm_extents, ONLY: extent, extent_shape, extent_size, is_contained_in, &
       extent_start
  USE iso_c_binding, ONLY: c_ptr, c_f_pointer
  USE mpi
  IMPLICIT NONE
  PRIVATE
  !INCLUDE 'ftype_size.inc'
  INTEGER, PARAMETER :: max_rank = 7
  INTEGER, PARAMETER :: extent_mp = MPI_2INTEGER
  INTEGER, PARAMETER :: ppm_address_kind = mpi_address_kind
  INTEGER(ppm_address_kind), PARAMETER :: ppm_maximum_alignment = 64

  !> describes one array-like global data structure to be distributed
  TYPE global_array_desc
    !> Fortran array rank of array
    INTEGER :: a_rank
    !> rect(1:a_rank) describes bounds of global array
    TYPE(extent) :: rect(max_rank)
    !> MPI datatype that describes the elements of an array
    INTEGER :: element_dt
  END TYPE global_array_desc

  INTEGER, PARAMETER :: ppm_default_comm = mpi_comm_world

  INTEGER, PARAMETER :: &
       not_exposed = 0, &
       exposed = 1

  INTEGER, PARAMETER, PUBLIC :: &
       !> in this mode calls to dist_mult_array_get will immediately
       !! retrieve the requested value
       sync_mode_passive_target = 0, &
       !> in this mode calls to dist_mult_array_get will result in
       !! the passed variable to become defined only after the next call
       !! to dist_mult_array_unexpose
       sync_mode_active_target = 1

  TYPE dm_array_cache_entry
    INTEGER(mpi_address_kind) :: base
    INTEGER :: composite_dt
    INTEGER :: rank
    INTEGER :: access_stamp
    INTEGER(mpi_address_kind), ALLOCATABLE :: offset(:)
  END TYPE dm_array_cache_entry

  TYPE dtype_info
    LOGICAL :: is_basic_type
    INTEGER(mpi_address_kind) :: extent
  END TYPE dtype_info

  !> dist_mult_array describes a global array where each rank only
  !! holds a part of the whole and this object can be used to access values
  !! stored on other ranks.
  TYPE dist_mult_array
    PRIVATE
    !> number of arrays that are distributed
    INTEGER :: num_sub_arrays
    !> Per distributed array information on global array shape and contents.\n
    !! The size of this array is 1:num_sub_arrays.
    TYPE(global_array_desc), ALLOCATABLE :: sub_arrays_global_desc(:)
    !> Stores local rectangular slices/chunks for each MPI rank.\n
    !! local_chunks(r, s, m) gives extent of rank r of sub-array s on MPI rank m
    TYPE(extent), ALLOCATABLE :: local_chunks(:, :, :)
    !> MPI communicator of all ranks holding the distributed array
    INTEGER :: comm
    !> size of above communicator
    INTEGER :: comm_size
    !> rank within above communicator
    INTEGER :: comm_rank
    !> MPI RDMA window associated with distributed array
    INTEGER :: win
    !> exposure status (either not_exposed or exposed)
    INTEGER :: exposure_status
    !> synchronization model (passive or active target)
    INTEGER :: sync_mode
    !> base extent of datatype of each sub-array, shape = (/ num_sub_arrays /)
    TYPE(dtype_info), ALLOCATABLE :: dt_info(:)
    !> stamp of most recent remote access
    INTEGER :: access_stamp, valid_stamp
    !> maximal window size over all members
    INTEGER(mpi_address_kind) :: max_win_size
    !> cache of remote data (index 0 contains local chunks)
    TYPE(dm_array_cache_entry), ALLOCATABLE :: cache(:)
  END TYPE dist_mult_array

#ifdef __SXorVH__
  INTEGER :: old_dist = 0
#endif
#ifdef __BLOCK_GET__
  TYPE iptr
    INTEGER(i4), ALLOCATABLE, DIMENSION(:) :: lst
    INTEGER :: maxelem
  END TYPE iptr 
  TYPE i8ptr
    INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: lst
  END TYPE i8ptr 
  TYPE dpptr
    REAL(dp), ALLOCATABLE, DIMENSION(:) :: lst
  END TYPE dpptr 
  INTEGER :: comm_0_cnt, comm_0_cnt_dp
  INTEGER,ALLOCATABLE :: comm_cnt(:), comm_cnt_dp(:)
  TYPE(iptr), ALLOCATABLE :: comm_lst(:)
  TYPE(iptr), ALLOCATABLE :: comm_lst_dp(:)
  TYPE(iptr),ALLOCATABLE :: recv_array(:)
  TYPE(dpptr),ALLOCATABLE :: recv_array_dp(:)
  INTEGER, ALLOCATABLE :: ones(:)
  INTEGER, ALLOCATABLE :: proc_lst(:,:,:), proc_lst_dp(:,:,:)
  INTEGER, ALLOCATABLE :: ind_lst(:,:,:), ind_lst_dp(:,:,:)
  INTEGER, ALLOCATABLE :: index_int_array(:)
  INTEGER, ALLOCATABLE :: index_dbl_array(:)
  ! variables for a more intelligent allocation process
  INTEGER :: lstblksizei4 = 100
  INTEGER :: lstblksizedp = 100
  
  PUBLIC :: recv_array, proc_lst, ind_lst
  PUBLIC :: recv_array_dp, proc_lst_dp, ind_lst_dp
#endif

  PUBLIC ppm_real_dp, ppm_real_sp, ppm_int, ppm_int_i8, ppm_bool

  PUBLIC :: dist_mult_array, global_array_desc
  PUBLIC :: dist_mult_array_new, dist_mult_array_delete
  PUBLIC :: dist_mult_array_local_ptr, dist_mult_array_get
  PUBLIC :: dist_mult_array_expose, dist_mult_array_unexpose
  PUBLIC :: dist_mult_array_rma_sync
#ifdef __BLOCK_GET__
  PUBLIC :: dist_mult_array_get_i4_blk, dist_mult_init_blk_comm, dist_mult_do_blk_comm, dist_mult_end_blk_comm
  PUBLIC :: dist_mult_array_get_dp_blk, dist_mult_init_blk_comm_dp, dist_mult_do_blk_comm_dp, dist_mult_end_blk_comm_dp
#endif

  !> @brief get POINTER to local portion of sub-array
  !!
  !! Obviously the rank of the POINTER argument has to match the dimensions
  !! used to specify the local part in the initialization of the
  !! distributed multi-array.
  !! In case the local part if of a non-intrinsic type, the user will
  !! need to get a \a TYPE(c_ptr) variable first and then convert to
  !! the corresponding type via \a C_F_POINTER.
  !!
  !! To change the data at positions i_l,j_l for a rank 2 sub-array of
  !! type INTEGER at index 5 use code along the following lines:
  !! @code
  !! TYPE(dist_mult_array) :: dma
  !! INTEGER, POINTER :: a5p(:, :)
  !! CALL dist_mult_array_local_ptr(dma, 5, a5p)
  !! CALL dist_mult_array_unexpose(dma)
  !! a5p(i_l, j_l) = i * j ! and other modifications of a5p...
  !! CALL dist_mult_array_expose(dma)
  !! @endcode
  INTERFACE dist_mult_array_local_ptr
    MODULE PROCEDURE dist_mult_array_local_ptr_c
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_1d
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_2d
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_3d
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_4d
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_5d
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_6d
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_7d
    MODULE PROCEDURE dist_mult_array_local_ptr_i8_1d
    MODULE PROCEDURE dist_mult_array_local_ptr_i8_2d
    MODULE PROCEDURE dist_mult_array_local_ptr_i8_3d
    MODULE PROCEDURE dist_mult_array_local_ptr_i8_4d
    MODULE PROCEDURE dist_mult_array_local_ptr_i8_5d
    MODULE PROCEDURE dist_mult_array_local_ptr_i8_6d
    MODULE PROCEDURE dist_mult_array_local_ptr_i8_7d
    MODULE PROCEDURE dist_mult_array_local_ptr_l_1d
    MODULE PROCEDURE dist_mult_array_local_ptr_l_2d
    MODULE PROCEDURE dist_mult_array_local_ptr_l_3d
    MODULE PROCEDURE dist_mult_array_local_ptr_l_4d
    MODULE PROCEDURE dist_mult_array_local_ptr_l_5d
    MODULE PROCEDURE dist_mult_array_local_ptr_l_6d
    MODULE PROCEDURE dist_mult_array_local_ptr_l_7d
    MODULE PROCEDURE dist_mult_array_local_ptr_sp_1d
    MODULE PROCEDURE dist_mult_array_local_ptr_sp_2d
    MODULE PROCEDURE dist_mult_array_local_ptr_sp_3d
    MODULE PROCEDURE dist_mult_array_local_ptr_sp_4d
    MODULE PROCEDURE dist_mult_array_local_ptr_sp_5d
    MODULE PROCEDURE dist_mult_array_local_ptr_sp_6d
    MODULE PROCEDURE dist_mult_array_local_ptr_sp_7d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_1d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_2d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_3d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_4d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_5d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_6d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_7d
  END INTERFACE dist_mult_array_local_ptr

  !> @brief Get value out of distributed multi-array, independent of
  !! rank the data resides on.
  !!
  !! For example, to query the data assigned in the example for
  !! dist_mult_array_local_ptr, use a query like this:
  !! @code
  !! TYPE(dist_mult_array) :: dma
  !! INTEGER :: i5
  !! CALL dist_mult_array_expose(dma)
  !! CALL dist_mult_array_get(dma, 5, (/ i_l, j_l /), i5)
  !! @endcode
  !! In case @a dma was initialized with sync_mode=sync_mode_active_target,
  !! the value of i5 must not be accessed before calling a synchronizing
  !! routine, either
  !! @code
  !! CALL dist_mult_array_unexpose(dma)
  !! @endcode
  !! or
  !! @code
  !! CALL dist_mult_array_rma_sync(dma)
  !! @endcode
  INTERFACE dist_mult_array_get
    MODULE PROCEDURE dist_mult_array_get_i4
    MODULE PROCEDURE dist_mult_array_get_i8
    MODULE PROCEDURE dist_mult_array_get_l
    MODULE PROCEDURE dist_mult_array_get_sp
    MODULE PROCEDURE dist_mult_array_get_dp
  END INTERFACE dist_mult_array_get

CONTAINS
  !> @brief create distributed multi-array data structure
  !!
  !! The resulting data type represents a number of arrays distributed
  !! over the ranks of the communicator passed to this function.
  !! @return initialized dist_mult_array structure
  !! @remark This operation is collective for all MPI ranks in @a comm
  FUNCTION dist_mult_array_new(sub_arrays, local_chunk, comm, &
       cache_size, sync_mode) &
       RESULT(dm_array)
    !> gives number (by its size), array ranks, data types
    !! and bounds of each distributed global array. The bounds are represented
    !! by an @ref extent type that stores start and size.
    TYPE(global_array_desc), INTENT(in) :: sub_arrays(:)
    !> local_chunk(i, j) describes for dimension i
    !! of sub_array j of the global arrays the local part available on
    !! this MPI rank.\n Only contiguous local parts are possible.\n
    !! <tt>SIZE(local_chunk, 2)</tt> must match
    !! <tt>SIZE(sub_arrays)</tt> and
    !! <tt>SHAPE(local_chunk) = (/ max_rank or more, num_sub_arrays /)</tt>
    TYPE(extent), INTENT(in) :: local_chunk(:, :)
    !> MPI communicator for which this data structure is
    !! collectively created
    INTEGER, INTENT(in) :: comm
    !> number of ranks to cache remote local parts of
    INTEGER, OPTIONAL, INTENT(in) :: cache_size
    !> switch synchronization modes, either sync_mode_passive_target
    !! (the default) or sync_mode_active_target_mode.  For
    !! sync_mode=sync_mode_active_target_mode, execution of RMA is
    !! deferred until the next synchronizing call
    !! (dist_mult_array_unexpose or dist_mult_array_rma_sync).
    !! (see @a dist_mult_array_get for example)
    INTEGER, OPTIONAL, INTENT(in) :: sync_mode
    TYPE(dist_mult_array) :: dm_array

    INTEGER :: ierror, num_sub_arrays, max_sub_array_rank, &
         comm_rank, comm_size, i, rank, dtype_prev
    INTEGER(mpi_address_kind) :: remote_win_size, max_remote_win_size
    INTEGER :: cache_size_, sync_mode_, win_create_info

    IF (PRESENT(sync_mode)) THEN
      sync_mode_ = sync_mode
      CALL assertion(sync_mode == IAND(sync_mode, 1), &
           __FILE__, &
           __LINE__, "sync_mode must be one of sync_mode_passive_target or &
           &sync_mode_active_target!")
    ELSE
      sync_mode_ = sync_mode_passive_target
    END IF


    num_sub_arrays = SIZE(sub_arrays)
    dm_array%num_sub_arrays = num_sub_arrays
    max_sub_array_rank = MAXVAL(sub_arrays%a_rank)
    CALL mpi_comm_size(comm, comm_size, ierror)
    dm_array%comm_size = comm_size
    CALL handle_mpi_error(ierror, comm, __LINE__)
    ALLOCATE(dm_array%sub_arrays_global_desc(num_sub_arrays), &
         dm_array%local_chunks(max_sub_array_rank, num_sub_arrays, &
         0:comm_size-1), &
         dm_array%dt_info(num_sub_arrays))
    dm_array%sub_arrays_global_desc = sub_arrays
    CALL mpi_comm_rank(comm, comm_rank, ierror)
    CALL handle_mpi_error(ierror, comm, __LINE__)
    dm_array%comm_rank = comm_rank
    ! copy local chunk bounds to all ranks
    dm_array%local_chunks(1:max_sub_array_rank, :, comm_rank) &
         = local_chunk(1:max_sub_array_rank, :)
    CALL mpi_allgather(mpi_in_place, 0, mpi_datatype_null, &
         dm_array%local_chunks, max_sub_array_rank * num_sub_arrays, &
         extent_mp, comm, ierror)
    CALL handle_mpi_error(ierror, comm, __LINE__)
    CALL mpi_comm_dup(comm, dm_array%comm, ierror)
    CALL handle_mpi_error(ierror, comm, __LINE__)

    dm_array%sync_mode = sync_mode_
    IF (sync_mode_ == sync_mode_passive_target) THEN
      IF (PRESENT(cache_size)) THEN
        cache_size_ = cache_size
      ELSE
        ! todo: this heuristic needs to improved!
        cache_size_ = MIN(CEILING(SQRT(REAL(comm_size))), max_sub_array_rank)
      END IF
      ALLOCATE(dm_array%cache(0:cache_size_))
      ! stamps are used to establish data age
      dm_array%access_stamp = 0
      dm_array%valid_stamp = 0
      ! initialize cache
      DO i = 0, cache_size_
        dm_array%cache(i)%rank = -1
        dm_array%cache(i)%base = 0
        dm_array%cache(i)%composite_dt = mpi_datatype_null
        dm_array%cache(i)%access_stamp = -1
        ALLOCATE(dm_array%cache(i)%offset(num_sub_arrays))
      END DO
      win_create_info = mpi_info_null
    ELSE ! sync_mode == sync_mode_active_target
      CALL mpi_info_create(win_create_info, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
      CALL mpi_info_set(win_create_info, 'no_locks', 'true', ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
      ALLOCATE(dm_array%cache(0:0))
      ALLOCATE(dm_array%cache(0)%offset(num_sub_arrays))
    END IF

    CALL get_extra_dtype_info
    CALL init_local_mem

    IF (sync_mode_ == sync_mode_passive_target) THEN
      ! establish maximal memory size of remote chunks to use in allocation of
      ! memory for caches
      max_remote_win_size = 0_mpi_address_kind
      DO rank = 0, comm_size - 1
        remote_win_size = 0_mpi_address_kind
        dtype_prev = mpi_datatype_null
        DO i = 1, num_sub_arrays
          remote_win_size = align_addr(remote_win_size, &
               &  dtype_prev, dm_array%sub_arrays_global_desc(i)%element_dt, &
               &  dm_array%dt_info(i)) &
               + dm_array%dt_info(i)%extent &
               * INT(extent_size(dm_array%local_chunks(1: &
               &                      sub_arrays(i)%a_rank, i, rank)), &
               &     mpi_address_kind)
          dtype_prev = dm_array%sub_arrays_global_desc(i)%element_dt
        END DO
        IF (remote_win_size > max_remote_win_size) &
             max_remote_win_size = remote_win_size
      END DO
      dm_array%max_win_size = max_remote_win_size
      CALL mpi_win_lock(mpi_lock_exclusive, dm_array%comm_rank, &
           mpi_mode_nocheck, dm_array%win, ierror)
    ELSE
      CALL mpi_info_free(win_create_info, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    END IF
  CONTAINS
    !> get extent for each sub_array element datatype
    SUBROUTINE get_extra_dtype_info
      INTEGER :: i
      INTEGER(mpi_address_kind) :: lb

      DO i = 1, num_sub_arrays
        CALL mpi_type_get_extent(sub_arrays(i)%element_dt, lb, &
             dm_array%dt_info(i)%extent, ierror)
        CALL handle_mpi_error(ierror, comm, __LINE__)
        dm_array%dt_info(i)%is_basic_type &
             = is_basic_type(sub_arrays(i)%element_dt, comm)
      END DO
    END SUBROUTINE get_extra_dtype_info

    !> create memory to store local parts in
    SUBROUTINE init_local_mem
      INTEGER(mpi_address_kind) :: local_win_size
      TYPE(c_ptr) :: baseptr_c
      INTEGER, POINTER :: baseptr

      CALL compute_cache_addr(dm_array%cache(0), sub_arrays, local_chunk, &
           dm_array%dt_info, local_win_size)
      CALL mpi_alloc_mem(local_win_size, mpi_info_null, &
           dm_array%cache(0)%base, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)

      baseptr_c = TRANSFER(dm_array%cache(0)%base, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
      CALL mpi_win_create(baseptr, local_win_size, 1, win_create_info, &
         dm_array%comm, dm_array%win, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
      dm_array%exposure_status = not_exposed
      dm_array%cache(0)%rank = comm_rank
    END SUBROUTINE init_local_mem
  END FUNCTION dist_mult_array_new

  !> @brief destruct dist_mult_array data type
  !! @remark since an MPI window is freed, this routine is collective for
  !! all processes which participated in the corresponding call to
  !! @ref dist_mult_array_new
  !! @param[inout] dm_array distributed multi-array to delete
  SUBROUTINE dist_mult_array_delete(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array

    INTEGER :: ierror, i

    IF (dm_array%exposure_status == exposed) THEN
      CALL dist_mult_array_unexpose(dm_array)
    END IF
    IF (dm_array%sync_mode == sync_mode_passive_target) &
         CALL mpi_win_unlock(dm_array%comm_rank, dm_array%win, ierror)
    CALL mpi_win_free(dm_array%win, ierror)
    CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)

    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      DO i = LBOUND(dm_array%cache, 1), UBOUND(dm_array%cache, 1)
        CALL delete_cache_entry(dm_array%cache(i), dm_array%comm)
      END DO
    END IF
    CALL mpi_comm_free(dm_array%comm, ierror)
    CALL handle_mpi_error(ierror, ppm_default_comm, __LINE__)
    dm_array%comm = mpi_comm_null

    DEALLOCATE(dm_array%sub_arrays_global_desc, dm_array%local_chunks, &
         dm_array%dt_info)

  END SUBROUTINE dist_mult_array_delete

  !> free allocated memory of data cached for a single rank
  SUBROUTINE delete_cache_entry(cache_entry, comm)
    TYPE(dm_array_cache_entry), INTENT(inout) :: cache_entry
    INTEGER, INTENT(in) :: comm

    CALL free_cache_entry_dt(cache_entry, comm)
    IF (cache_entry%base /= 0) THEN
      CALL free_mpi_allocation(cache_entry%base, comm)
    END IF
  END SUBROUTINE delete_cache_entry

  !> free MPI datatype used to transfer remote rank memory
  SUBROUTINE free_cache_entry_dt(cache_entry, comm)
    TYPE(dm_array_cache_entry), INTENT(inout) :: cache_entry
    INTEGER, INTENT(in) :: comm

    INTEGER :: ierror

    IF (cache_entry%rank >= 0) THEN
      IF (cache_entry%composite_dt /= mpi_datatype_null) THEN
        CALL mpi_type_free(cache_entry%composite_dt, ierror)
        CALL handle_mpi_error(ierror, comm, __LINE__)
      END IF
    END IF
  END SUBROUTINE free_cache_entry_dt

  !> converts base address into POINTER that can be passed to mpi_free_mem
  SUBROUTINE free_mpi_allocation(address, comm)
    INTEGER(mpi_address_kind), INTENT(inout) :: address
    INTEGER, INTENT(in) :: comm

    TYPE(c_ptr) :: baseptr_c
    INTEGER, POINTER :: baseptr
    INTEGER :: ierror

    baseptr_c = TRANSFER(address, baseptr_c)
    CALL C_F_POINTER(baseptr_c, baseptr)
    CALL mpi_free_mem(baseptr, ierror)
    CALL handle_mpi_error(ierror, comm, __LINE__)
  END SUBROUTINE free_mpi_allocation

  !> computes offset each sub-array is stored at in cache memory
  !! optionally computes total size of corresponding window area
  SUBROUTINE compute_cache_addr(cache_entry, sub_arrays, chunk, dt_info, &
       win_size)
    TYPE(dm_array_cache_entry), INTENT(inout) :: cache_entry
    TYPE(global_array_desc), INTENT(in) :: sub_arrays(:)
    TYPE(extent), INTENT(in) :: chunk(:, :)
    TYPE(dtype_info), INTENT(in) :: dt_info(:)
    INTEGER(mpi_address_kind), OPTIONAL, INTENT(out) :: win_size

    INTEGER :: i, num_sub_arrays, dtprev
    INTEGER(mpi_address_kind) :: win_size_

    num_sub_arrays = SIZE(sub_arrays)
    win_size_ = 0_mpi_address_kind
    dtprev = mpi_datatype_null
    DO i = 1, num_sub_arrays
      win_size_ = align_addr(win_size_, dtprev, &
           sub_arrays(i)%element_dt, dt_info(i))
      cache_entry%offset(i) = win_size_
      win_size_ = win_size_ + dt_info(i)%extent &
           * INT(extent_size(chunk(1:sub_arrays(i)%a_rank, i)), &
           &     mpi_address_kind)
      dtprev = sub_arrays(i)%element_dt
    END DO
    IF (PRESENT(win_size)) win_size = win_size_
  END SUBROUTINE compute_cache_addr

  FUNCTION align_addr(addr, dtype_prev, dtype2align_for, &
       dt_info) RESULT(aligned_addr)
    INTEGER(mpi_address_kind), INTENT(in) :: addr
    INTEGER, INTENT(in) :: dtype_prev, dtype2align_for
    TYPE(dtype_info), INTENT(in) :: dt_info
    INTEGER(mpi_address_kind) :: aligned_addr

    IF (dtype_prev == mpi_datatype_null .OR. dtype_prev == dtype2align_for) THEN
      aligned_addr = addr
    ELSE IF (dt_info%is_basic_type) THEN
      ! align basic type to multiple of its extent
      aligned_addr = ((addr + dt_info%extent - 1) &
           / dt_info%extent) * dt_info%extent
    ELSE
      ! align non-basic type not at start of composite
      ! round up to next multiple of ppm_maximum_alignment
      IF (IAND(ppm_maximum_alignment, &
           ppm_maximum_alignment - 1_ppm_address_kind) &
           == 0_ppm_address_kind) THEN
        aligned_addr = IAND(addr + ppm_maximum_alignment - 1, &
             -ppm_maximum_alignment)
      ELSE
        aligned_addr = ((addr + ppm_maximum_alignment - 1) &
             / ppm_maximum_alignment) * ppm_maximum_alignment
      END IF
    END IF
  END FUNCTION align_addr

  !> determine if type is a basic datatype
  RECURSIVE FUNCTION is_basic_type(dtype, comm) RESULT(is_basic_type_)
    INTEGER, INTENT(in) :: dtype
    INTEGER :: comm
    LOGICAL :: is_basic_type_
    INTEGER(mpi_address_kind) :: cont_a(1)
    INTEGER :: num_int, num_addr, num_dt, combiner, ierror, &
         cont_i(1), cont_dt(1)
    CALL mpi_type_get_envelope(dtype, num_int, num_addr, num_dt, combiner, &
         ierror)
    CALL handle_mpi_error(ierror, comm, __LINE__)
    IF (combiner == mpi_combiner_named) THEN
      is_basic_type_ = .TRUE.
    ELSE IF (combiner == mpi_combiner_dup) THEN
      CALL mpi_type_get_contents(dtype, num_int, num_addr, num_dt, &
           cont_i, cont_a, cont_dt, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
      is_basic_type_ = is_basic_type(cont_dt(1), comm)
      CALL mpi_type_get_envelope(cont_dt(1), num_int, num_addr, num_dt, &
           combiner, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
      IF (.NOT. combiner == mpi_combiner_named) THEN
        CALL mpi_type_free(cont_dt(1), ierror)
        CALL handle_mpi_error(ierror, comm, __LINE__)
      END IF
    ELSE
      is_basic_type_ = .FALSE.
    END IF
  END FUNCTION is_basic_type

  !> retrieve TYPE(c_ptr) corresponding to a cached sub-array chunk
  SUBROUTINE dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    TYPE(c_ptr), INTENT(out) :: sub_array_ptr
    sub_array_ptr = TRANSFER(dm_array%cache(cache_idx)%base &
         + dm_array%cache(cache_idx)%offset(sub_array_idx), sub_array_ptr)
  END SUBROUTINE dist_mult_array_cache_ptr_c

  !> retrieve TYPE(c_ptr) to local chunk
  SUBROUTINE dist_mult_array_local_ptr_c(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    TYPE(c_ptr), INTENT(out) :: sub_array_ptr
    sub_array_ptr = TRANSFER(dm_array%cache(0)%base &
         + dm_array%cache(0)%offset(sub_array_idx), sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_c

  !> @brief make local data available to other ranks
  !!
  !! This call is collective for all ranks in the communicator used
  !! to create @a dm_array
  !!
  !! note: this enters an exposure epoch on the data shared via RMA
  !! local data must not be changed while the array is in
  !! exposed state
  SUBROUTINE dist_mult_array_expose(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array

    INTEGER, PARAMETER :: fence_assert = IOR(mpi_mode_noput, mpi_mode_noprecede)
    INTEGER :: ierror
    CALL assertion(dm_array%exposure_status == not_exposed, &
         __FILE__, &
         __LINE__, "distributed multi-array must not already be exposed!")
    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      CALL mpi_win_unlock(dm_array%comm_rank, dm_array%win, ierror)
      CALL mpi_barrier(dm_array%comm, ierror)
    ELSE ! dm_array%sync_mode == sync_mode_active_target
      CALL mpi_win_fence(fence_assert, dm_array%win, ierror)
    END IF
    CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
    dm_array%exposure_status = exposed
    dm_array%valid_stamp = dm_array%access_stamp + 1
  END SUBROUTINE dist_mult_array_expose

  !> @brief wait for all ranks to finish queries in current exposure epoch
  !!
  !! This call is collective for all ranks in the communicator used
  !! to create @a dm_array
  !!
  !! note: local data can only be changed after the distributed array
  !! was initially created or dist_mult_array_unexpose has been called
  SUBROUTINE dist_mult_array_unexpose(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array

    INTEGER :: fence_assert = IOR(IOR(mpi_mode_nostore, mpi_mode_noput), &
         mpi_mode_nosucceed)
    INTEGER :: ierror
    CALL assertion(dm_array%exposure_status == exposed, &
         __FILE__, &
         __LINE__, &
         "distributed multi-array must have previously been made accessible!")
    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      CALL mpi_barrier(dm_array%comm, ierror)
      CALL mpi_win_lock(mpi_lock_exclusive, dm_array%comm_rank, &
           mpi_mode_nocheck, dm_array%win, ierror)
    ELSE ! dm_array%sync_mode == sync_mode_active_target
      CALL mpi_win_fence(fence_assert, dm_array%win, ierror)
    END IF
    CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
    dm_array%exposure_status = not_exposed
  END SUBROUTINE dist_mult_array_unexpose

  !> @brief synchronize RMA updates only, ignore local updates
  !!
  !! This call is collective for all ranks in the communicator used
  !! to create @a dm_array.
  !! Also @a dm_array must be in exposed state.
  SUBROUTINE dist_mult_array_rma_sync(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array

    INTEGER :: fence_assert = IOR(mpi_mode_nostore, mpi_mode_noput)
    INTEGER :: ierror
    CALL assertion(dm_array%exposure_status == exposed, &
         __FILE__, &
         __LINE__, &
         "distributed multi-array must have previously been made accessible!")
    CALL assertion(dm_array%sync_mode == sync_mode_active_target, &
         __FILE__, &
         __LINE__, &
         "incorrect synchronization call")
    ! dm_array%sync_mode == sync_mode_active_target
    CALL mpi_win_fence(fence_assert, dm_array%win, ierror)
    CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
  END SUBROUTINE dist_mult_array_rma_sync

  !> map remote memory of rank to cache entry and return cache entry index
  FUNCTION dist_mult_array_cache_rank(dm_array, rank) RESULT(cache_entry)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: rank
    INTEGER :: cache_entry

    INTEGER :: i, oldest, stamp, oldest_stamp, ierror
    INTEGER :: component_count(dm_array%num_sub_arrays)
    INTEGER :: a_rank, num_sub_arrays

    IF (rank == dm_array%comm_rank) THEN
      cache_entry = 0
      RETURN
    END IF
    CALL assertion(dm_array%exposure_status == exposed, &
         __FILE__, &
         __LINE__, &
         "distributed multi-array must have previously been made accessible!")
    stamp = dm_array%access_stamp + 1
    oldest_stamp = stamp
    DO i = 1, UBOUND(dm_array%cache, 1)
      IF (dm_array%cache(i)%rank == rank) THEN
        IF (dm_array%cache(i)%access_stamp < dm_array%valid_stamp) &
             CALL cache_get(i, rank)
        dm_array%cache(i)%access_stamp = stamp
        dm_array%access_stamp = stamp
        cache_entry = i
        RETURN
      ELSE IF (dm_array%cache(i)%access_stamp < oldest_stamp) THEN
        oldest_stamp = dm_array%cache(i)%access_stamp
        oldest = i
      END IF
    END DO
    IF (dm_array%cache(oldest)%rank >= 0) THEN
      ! free cache resources here
      CALL free_cache_entry_dt(dm_array%cache(oldest), dm_array%comm)
    END IF
    ! populate free cache entry
    IF (dm_array%cache(oldest)%base == 0) THEN
      CALL mpi_alloc_mem(dm_array%max_win_size, mpi_info_null, &
           dm_array%cache(oldest)%base, ierror)
      CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
    END IF
    ! create datatype for mpi_get
    CALL compute_cache_addr(dm_array%cache(oldest), &
         dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, rank), dm_array%dt_info)

    num_sub_arrays = dm_array%num_sub_arrays
    DO i = 1, num_sub_arrays
      a_rank = dm_array%sub_arrays_global_desc(i)%a_rank
      component_count(i) = extent_size(dm_array%local_chunks(1:a_rank, i, rank))
    END DO
    CALL mpi_type_create_struct(num_sub_arrays, component_count, &
         dm_array%cache(oldest)%offset, &
         dm_array%sub_arrays_global_desc(:)%element_dt, &
         dm_array%cache(oldest)%composite_dt, ierror)
    CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
    CALL mpi_type_commit(dm_array%cache(oldest)%composite_dt, ierror)
    CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
    CALL cache_get(oldest, rank)

    dm_array%cache(oldest)%rank = rank

    dm_array%cache(oldest)%access_stamp = stamp
    dm_array%access_stamp = stamp

    cache_entry = oldest
  CONTAINS
    !> perform RMA get operation
    SUBROUTINE cache_get(cache_entry, rank)
      INTEGER, INTENT(in) :: cache_entry, rank

      INTEGER :: ierror
      TYPE(c_ptr) :: baseptr_c
      INTEGER, POINTER :: baseptr

      baseptr_c = TRANSFER(dm_array%cache(cache_entry)%base, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
      CALL mpi_win_lock(mpi_lock_shared, rank, mpi_mode_nocheck, &
           dm_array%win, ierror)
      CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
      CALL mpi_get(baseptr, 1, dm_array%cache(cache_entry)%composite_dt, &
           rank, 0_mpi_address_kind, 1,&
           dm_array%cache(cache_entry)%composite_dt, dm_array%win, ierror)
      CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
      CALL mpi_win_unlock(rank, dm_array%win, ierror)
      CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
    END SUBROUTINE cache_get
  END FUNCTION dist_mult_array_cache_rank

  FUNCTION dist_mult_array_coord2rank(dm_array, sub_array, coord) &
       RESULT(rank)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER :: rank, ref_rank

    INTEGER :: comm_rank, comm_size, max_dist, dist
#ifdef __SXorVH__
    INTEGER :: start_dist

    ! load dist that was successful in last last call
    start_dist = old_dist 
#endif

    comm_rank = dm_array%comm_rank
    comm_size = dm_array%comm_size
    max_dist = (comm_size + 1)/2
    ref_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
    DO dist = 0, max_dist
      IF (is_contained_in(coord, dm_array%local_chunks(1:ref_rank, sub_array, &
#ifndef __SXorVH__
           MOD(comm_rank + dist, comm_size)))) THEN
        rank = MOD(comm_rank + dist, comm_size)
#else
           MODULO(comm_rank + start_dist + dist, comm_size)))) THEN
        old_dist = start_dist + dist
        rank = MODULO(comm_rank + old_dist, comm_size)
#endif
        RETURN
      ELSE IF (is_contained_in(coord, dm_array%local_chunks(1:ref_rank, sub_array, &
#ifndef __SXorVH__
           MODULO(comm_rank - dist, comm_size)))) THEN
        rank = MODULO(comm_rank - dist, comm_size)
#else
           MODULO(comm_rank + start_dist - dist, comm_size)))) THEN
        old_dist = start_dist - dist
        rank = MODULO(comm_rank + old_dist, comm_size)
#endif
        RETURN
      END IF
    END DO
#ifdef __SXorVH__
! For debugging: second try without search optimization if first try was unsuccessful
    write(0,*) 'PROBLEM: repeated search in SR dist_mult_array_coord2rank because the first one was unsuccessful', &
     'comm_rank = ',comm_rank, 'comm_size = ', comm_size, 'max_dist = ', max_dist, 'ref_rank = ', ref_rank, 'p_pe =', p_pe
    DO dist = 0, max_dist
      IF (is_contained_in(coord, dm_array%local_chunks(1:ref_rank, sub_array, &
           MOD(comm_rank + dist, comm_size)))) THEN
        rank = MOD(comm_rank + dist, comm_size)
        RETURN
      ELSE IF (is_contained_in(coord, dm_array%local_chunks(1:ref_rank, sub_array, &
           MODULO(comm_rank - dist, comm_size)))) THEN
        rank = MODULO(comm_rank - dist, comm_size)
        RETURN
      END IF
    END DO
#endif
    rank = -1
    CALL abort_ppm("invalid coordinate", &
         __FILE__, &
         __LINE__, dm_array%comm)

  END FUNCTION dist_mult_array_coord2rank

  !> establish remote MPI rank holding data of sub-array at coordinate
  !> and fetch data if not already in-cache
  FUNCTION dist_mult_array_get_cache_idx(dm_array, sub_array, coord, &
       ref_extent) RESULT(cache_idx)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(mpi_address_kind), INTENT(in) :: ref_extent

    INTEGER :: ref_rank
    INTEGER :: cache_idx, ierror
    INTEGER(mpi_address_kind) :: lb, extent

    CALL assertion(sub_array >= 1 &
         .AND. sub_array <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    ref_rank = SIZE(coord)
    CALL assertion(ref_rank &
         == dm_array%sub_arrays_global_desc(sub_array)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion(is_contained_in(coord, &
         dm_array%sub_arrays_global_desc(sub_array)%rect(1:ref_rank)), &
         __FILE__, &
         __LINE__, "invalid coordinate")
    CALL mpi_type_get_extent(&
         dm_array%sub_arrays_global_desc(sub_array)%element_dt, &
         lb, extent, ierror)
    CALL handle_mpi_error(ierror, dm_array%comm, __LINE__)
    CALL assertion(extent == ref_extent, &
         __FILE__, &
         __LINE__, "extent mismatch")
    cache_idx = dist_mult_array_cache_rank(dm_array, &
         dist_mult_array_coord2rank(dm_array, sub_array, coord))
  END FUNCTION dist_mult_array_get_cache_idx



  ! see @ref dist_mult_array_get
  SUBROUTINE dist_mult_array_cache_val_i4(dm_array, sub_array_idx, &
       cache_idx, coord, v)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx, cache_idx
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i4), INTENT(out) :: v

    INTEGER(i4), POINTER :: p1(:), &
         &         p2(:,:), &
         &         p3(:,:,:), &
         &         p4(:,:,:,:), &
         &         p5(:,:,:,:,:), &
         &         p6(:,:,:,:,:,:), &
         &         p7(:,:,:,:,:,:,:)

    SELECT CASE (SIZE(coord))

    CASE (1)
      CALL dist_mult_array_cache_ptr_i4_1d(dm_array, sub_array_idx, &
           cache_idx, p1)
      v = p1(coord(1))
    CASE (2)
      CALL dist_mult_array_cache_ptr_i4_2d(dm_array, sub_array_idx, &
           cache_idx, p2)
      v = p2(coord(1),coord(2))
    CASE (3)
      CALL dist_mult_array_cache_ptr_i4_3d(dm_array, sub_array_idx, &
           cache_idx, p3)
      v = p3(coord(1),coord(2),coord(3))
    CASE (4)
      CALL dist_mult_array_cache_ptr_i4_4d(dm_array, sub_array_idx, &
           cache_idx, p4)
      v = p4(coord(1),coord(2),coord(3),coord(4))
    CASE (5)
      CALL dist_mult_array_cache_ptr_i4_5d(dm_array, sub_array_idx, &
           cache_idx, p5)
      v = p5(coord(1),coord(2),coord(3),coord(4),coord(5))
    CASE (6)
      CALL dist_mult_array_cache_ptr_i4_6d(dm_array, sub_array_idx, &
           cache_idx, p6)
      v = p6(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6))
    CASE (7)
      CALL dist_mult_array_cache_ptr_i4_7d(dm_array, sub_array_idx, &
           cache_idx, p7)
      v = p7(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6),coord(7))
    CASE default
      CALL abort_ppm("invalid array rank", &
           __FILE__, &
           __LINE__, dm_array%comm)
    END SELECT
  END SUBROUTINE dist_mult_array_cache_val_i4


  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i4_1d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(1)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:1, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i4_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i4_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:)

    CALL dist_mult_array_cache_ptr_i4_1d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i4_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i4_2d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(2)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:2, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i4_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i4_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:)

    CALL dist_mult_array_cache_ptr_i4_2d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i4_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i4_3d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(3)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:3, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i4_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i4_3d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:)

    CALL dist_mult_array_cache_ptr_i4_3d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i4_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i4_4d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(4)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:4, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i4_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i4_4d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:)

    CALL dist_mult_array_cache_ptr_i4_4d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i4_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i4_5d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(5)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:5, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i4_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i4_5d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_i4_5d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i4_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i4_6d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(6)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:6, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i4_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i4_6d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_i4_6d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i4_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i4_7d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(7)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:7, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(7, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i4_7d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i4_7d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_i4_7d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i4_7d


  SUBROUTINE dist_mult_array_get_deferred_i4(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i4), INTENT(out) :: v

    INTEGER(mpi_address_kind) :: byte_offset, ofs_factor
    INTEGER :: coord_base(7), a_rank, i
    TYPE(dm_array_cache_entry) :: cache_entry
    INTEGER :: src_comm_rank, comm, dt, ierror
    TYPE(c_ptr) :: baseptr_c
    INTEGER(i4), POINTER :: baseptr

    comm = dm_array%comm
    src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord)
    ALLOCATE(cache_entry%offset(dm_array%num_sub_arrays))
    CALL compute_cache_addr(cache_entry, dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, src_comm_rank), dm_array%dt_info)
    a_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
    coord_base(1:a_rank) &
         = dm_array%local_chunks(1:a_rank, sub_array, src_comm_rank)%first
    byte_offset = cache_entry%offset(sub_array)
    ofs_factor = dm_array%dt_info(sub_array)%extent
    DO i = 1, a_rank
      byte_offset = byte_offset + (coord(i) - coord_base(i)) * ofs_factor
      ofs_factor = ofs_factor &
           * INT(dm_array%local_chunks(i, sub_array, src_comm_rank)%size, &
           &     mpi_address_kind)
    END DO
    dt = dm_array%sub_arrays_global_desc(sub_array)%element_dt

    IF (src_comm_rank == dm_array%comm_rank) THEN
      baseptr_c = TRANSFER(dm_array%cache(0)%base + byte_offset, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
      CALL mpi_sendrecv(baseptr, 1, dt, 0, 0, v, 1, dt, 0, 0, &
           mpi_comm_self, mpi_status_ignore, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    ELSE
      CALL mpi_get(v, 1, dt, src_comm_rank, byte_offset, 1, dt, dm_array%win, &
           ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    END IF

  END SUBROUTINE dist_mult_array_get_deferred_i4

#ifdef __BLOCK_GET__
  !NEC: replacement for dist_mult_array_get_deferred_i4: where all data is
  !stored instead of get directly
  SUBROUTINE dist_mult_array_get_i4_blk(dm_array, sub_array, coord, v, proc, ind)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i4), INTENT(out) :: v
    INTEGER, INTENT(inout) :: proc 
    INTEGER, INTENT(inout) :: ind 

    INTEGER(mpi_address_kind) :: byte_offset, ofs_factor
    INTEGER :: coord_base(7), a_rank, i
    TYPE(dm_array_cache_entry) :: cache_entry
    INTEGER :: src_comm_rank, comm, dt, ierror
    TYPE(c_ptr) :: baseptr_c
    INTEGER(i4), POINTER :: baseptr
    INTEGER :: omaxsize
    INTEGER, ALLOCATABLE, DIMENSION(:) :: tmparr
    comm = dm_array%comm
    src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord)
    ALLOCATE(cache_entry%offset(dm_array%num_sub_arrays))
    CALL compute_cache_addr(cache_entry, dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, src_comm_rank), dm_array%dt_info)
    a_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
!$NEC novector
    coord_base(1:a_rank) &
         = dm_array%local_chunks(1:a_rank, sub_array, src_comm_rank)%first
    byte_offset = cache_entry%offset(sub_array)
    ofs_factor = dm_array%dt_info(sub_array)%extent
!$NEC novector
    DO i = 1, a_rank
      byte_offset = byte_offset + (coord(i) - coord_base(i)) * ofs_factor
      ofs_factor = ofs_factor &
           * INT(dm_array%local_chunks(i, sub_array, src_comm_rank)%size, &
           &     mpi_address_kind)
    END DO
    dt = dm_array%sub_arrays_global_desc(sub_array)%element_dt

    IF (src_comm_rank == dm_array%comm_rank) THEN
      comm_0_cnt = comm_0_cnt + 1
      baseptr_c = TRANSFER(dm_array%cache(0)%base + byte_offset, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
#ifndef __COMM_OPT__
      CALL mpi_sendrecv(baseptr, 1, dt, 0, 0, v, 1, dt, 0, 0, &
           mpi_comm_self, mpi_status_ignore, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
#else
      v = baseptr
#endif
    ELSE
      comm_cnt(src_comm_rank) = comm_cnt(src_comm_rank) + 1
      IF (comm_cnt(src_comm_rank) > comm_lst(src_comm_rank)%maxelem) THEN 
        !NEC: reallocation procedure in order to save some memory; tuned, to be efficient
        IF (ALLOCATED(comm_lst(src_comm_rank)%lst)) THEN
          omaxsize = comm_lst(src_comm_rank)%maxelem
          ALLOCATE(tmparr(omaxsize))
          tmparr = comm_lst(src_comm_rank)%lst
          DEALLOCATE(comm_lst(src_comm_rank)%lst)
          comm_lst(src_comm_rank)%maxelem = INT(1.41*omaxsize)
          ALLOCATE(comm_lst(src_comm_rank)%lst(comm_lst(src_comm_rank)%maxelem))
          DO i = 1, omaxsize
            comm_lst(src_comm_rank)%lst(i) = tmparr(i)
          END DO
          DEALLOCATE(tmparr)
        ELSE
          comm_lst(src_comm_rank)%maxelem = lstblksizei4
          ALLOCATE(comm_lst(src_comm_rank)%lst(comm_lst(src_comm_rank)%maxelem))
        END IF
      END IF
     !NEC: store offset and return proc and the access number to this proc
      comm_lst(src_comm_rank)%lst(comm_cnt(src_comm_rank)) = byte_offset
      proc = src_comm_rank
      ind = comm_cnt(src_comm_rank)
    END IF

  END SUBROUTINE dist_mult_array_get_i4_blk

  !NEC: initializes data structure for blocked communication
  SUBROUTINE dist_mult_init_blk_comm(dim1, dim2, dim3, dm_array)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: dim1, dim2, dim3
    !comm_max = dim1 * dim2 * dim3 
    ALLOCATE(comm_cnt(0:dm_array%comm_size-1), comm_lst(0:dm_array%comm_size-1), &
             recv_array(0:dm_array%comm_size-1), index_int_array(0:dm_array%comm_size-1))
    comm_lst(:)%maxelem = 0
    comm_0_cnt = 0
    comm_cnt = 0
    ALLOCATE(proc_lst(dim1, dim2, dim3), ind_lst(dim1, dim2,dim3))
    proc_lst = -1
  END SUBROUTINE dist_mult_init_blk_comm

  !NEC: executes blocked communcation
  SUBROUTINE dist_mult_do_blk_comm(dm_array)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER :: i, j, ierr
    INTEGER, SAVE :: tsize = 0

   !NEC: Only allocate "ones" once to save time
    ALLOCATE(ones(MAXVAL(comm_cnt)))
    ones = 1
   !NEC: Store MPI type size assuming it is not changing on a process during a run
    IF (tsize == 0) CALL MPI_TYPE_SIZE(MPI_INTEGER, tsize, ierr)
    DO i = 0, dm_array%comm_size-1
      IF (ALLOCATED(comm_lst(i)%lst)) THEN
        !NEC: Make byte offset list into type offset list to be portable
        DO j = 1,comm_cnt(i)
          comm_lst(i)%lst(j) = comm_lst(i)%lst(j) / tsize
        END DO
        call MPI_TYPE_INDEXED(comm_cnt(i), ones, comm_lst(i)%lst, MPI_INTEGER, index_int_array(i), ierr)
        call MPI_TYPE_COMMIT(index_int_array(i), ierr)
        ALLOCATE(recv_array(i)%lst(comm_cnt(i)))
        call MPI_GET(recv_array(i)%lst, comm_cnt(i), MPI_INTEGER, &
                     i, 0_mpi_address_kind, 1, index_int_array(i), &
                     dm_array%win, ierr) 
      END IF
    END DO
    DEALLOCATE(ones)
  END SUBROUTINE dist_mult_do_blk_comm

  !NEC: frees data structure for blocked communication
  SUBROUTINE dist_mult_end_blk_comm(dm_array)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER :: i, ierr

    DO i = 0, dm_array%comm_size-1
      IF (ALLOCATED(comm_lst(i)%lst)) THEN
        DEALLOCATE(comm_lst(i)%lst)
        DEALLOCATE(recv_array(i)%lst)
        call MPI_TYPE_FREE(index_int_array(i), ierr)
      END IF
    END DO
    !NEC: Make the next selected array size the largest needed, to save reallocations, thus time
    lstblksizei4 = MAX(MAXVAL(comm_cnt),100)
    DEALLOCATE(comm_cnt, comm_lst, recv_array, index_int_array)
    DEALLOCATE(proc_lst, ind_lst)
  END SUBROUTINE dist_mult_end_blk_comm

#endif

  SUBROUTINE dist_mult_array_get_i4(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i4), INTENT(out) :: v

    INTEGER :: cache_idx

    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      cache_idx = dist_mult_array_get_cache_idx(dm_array, sub_array, coord, &
           INT(mp_i4_extent, mpi_address_kind))
      CALL dist_mult_array_cache_val_i4(dm_array, sub_array, &
           cache_idx, coord, v)
    ELSE ! dm_array%sync_mode == sync_mode_active_target
      CALL dist_mult_array_get_deferred_i4(dm_array, sub_array, coord, v)
    END IF
  END SUBROUTINE dist_mult_array_get_i4


  ! see @ref dist_mult_array_get
  SUBROUTINE dist_mult_array_cache_val_i8(dm_array, sub_array_idx, &
       cache_idx, coord, v)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx, cache_idx
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i8), INTENT(out) :: v

    INTEGER(i8), POINTER :: p1(:), &
         &         p2(:,:), &
         &         p3(:,:,:), &
         &         p4(:,:,:,:), &
         &         p5(:,:,:,:,:), &
         &         p6(:,:,:,:,:,:), &
         &         p7(:,:,:,:,:,:,:)

    SELECT CASE (SIZE(coord))

    CASE (1)
      CALL dist_mult_array_cache_ptr_i8_1d(dm_array, sub_array_idx, &
           cache_idx, p1)
      v = p1(coord(1))
    CASE (2)
      CALL dist_mult_array_cache_ptr_i8_2d(dm_array, sub_array_idx, &
           cache_idx, p2)
      v = p2(coord(1),coord(2))
    CASE (3)
      CALL dist_mult_array_cache_ptr_i8_3d(dm_array, sub_array_idx, &
           cache_idx, p3)
      v = p3(coord(1),coord(2),coord(3))
    CASE (4)
      CALL dist_mult_array_cache_ptr_i8_4d(dm_array, sub_array_idx, &
           cache_idx, p4)
      v = p4(coord(1),coord(2),coord(3),coord(4))
    CASE (5)
      CALL dist_mult_array_cache_ptr_i8_5d(dm_array, sub_array_idx, &
           cache_idx, p5)
      v = p5(coord(1),coord(2),coord(3),coord(4),coord(5))
    CASE (6)
      CALL dist_mult_array_cache_ptr_i8_6d(dm_array, sub_array_idx, &
           cache_idx, p6)
      v = p6(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6))
    CASE (7)
      CALL dist_mult_array_cache_ptr_i8_7d(dm_array, sub_array_idx, &
           cache_idx, p7)
      v = p7(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6),coord(7))
    CASE default
      CALL abort_ppm("invalid array rank", &
           __FILE__, &
           __LINE__, dm_array%comm)
    END SELECT
  END SUBROUTINE dist_mult_array_cache_val_i8


  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i8_1d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(1)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:1, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i8_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i8_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:)

    CALL dist_mult_array_cache_ptr_i8_1d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i8_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i8_2d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(2)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:2, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i8_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i8_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:)

    CALL dist_mult_array_cache_ptr_i8_2d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i8_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i8_3d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(3)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:3, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i8_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i8_3d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:)

    CALL dist_mult_array_cache_ptr_i8_3d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i8_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i8_4d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(4)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:4, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i8_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i8_4d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:)

    CALL dist_mult_array_cache_ptr_i8_4d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i8_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i8_5d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(5)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:5, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i8_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i8_5d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_i8_5d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i8_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i8_6d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(6)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:6, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i8_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i8_6d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_i8_6d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i8_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_i8_7d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(7)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:7, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(7, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_i8_7d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_i8_7d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i8), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_i8_7d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_i8_7d


  SUBROUTINE dist_mult_array_get_deferred_i8(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i8), INTENT(out) :: v

    INTEGER(mpi_address_kind) :: byte_offset, ofs_factor
    INTEGER :: coord_base(7), a_rank, i
    TYPE(dm_array_cache_entry) :: cache_entry
    INTEGER :: src_comm_rank, comm, dt, ierror
    TYPE(c_ptr) :: baseptr_c
    INTEGER(i8), POINTER :: baseptr

    comm = dm_array%comm
    src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord)
    ALLOCATE(cache_entry%offset(dm_array%num_sub_arrays))
    CALL compute_cache_addr(cache_entry, dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, src_comm_rank), dm_array%dt_info)
    a_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
    coord_base(1:a_rank) &
         = dm_array%local_chunks(1:a_rank, sub_array, src_comm_rank)%first
    byte_offset = cache_entry%offset(sub_array)
    ofs_factor = dm_array%dt_info(sub_array)%extent
    DO i = 1, a_rank
      byte_offset = byte_offset + (coord(i) - coord_base(i)) * ofs_factor
      ofs_factor = ofs_factor &
           * INT(dm_array%local_chunks(i, sub_array, src_comm_rank)%size, &
           &     mpi_address_kind)
    END DO
    dt = dm_array%sub_arrays_global_desc(sub_array)%element_dt

    IF (src_comm_rank == dm_array%comm_rank) THEN
      baseptr_c = TRANSFER(dm_array%cache(0)%base + byte_offset, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
      CALL mpi_sendrecv(baseptr, 1, dt, 0, 0, v, 1, dt, 0, 0, &
           mpi_comm_self, mpi_status_ignore, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    ELSE
      CALL mpi_get(v, 1, dt, src_comm_rank, byte_offset, 1, dt, dm_array%win, &
           ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    END IF

  END SUBROUTINE dist_mult_array_get_deferred_i8

  SUBROUTINE dist_mult_array_get_i8(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i8), INTENT(out) :: v

    INTEGER :: cache_idx

    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      cache_idx = dist_mult_array_get_cache_idx(dm_array, sub_array, coord, &
           INT(mp_i8_extent, mpi_address_kind))
      CALL dist_mult_array_cache_val_i8(dm_array, sub_array, &
           cache_idx, coord, v)
    ELSE ! dm_array%sync_mode == sync_mode_active_target
      CALL dist_mult_array_get_deferred_i8(dm_array, sub_array, coord, v)
    END IF
  END SUBROUTINE dist_mult_array_get_i8


  ! see @ref dist_mult_array_get
  SUBROUTINE dist_mult_array_cache_val_l(dm_array, sub_array_idx, &
       cache_idx, coord, v)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx, cache_idx
    INTEGER, INTENT(in) :: coord(:)
    LOGICAL, INTENT(out) :: v

    LOGICAL, POINTER :: p1(:), &
         &         p2(:,:), &
         &         p3(:,:,:), &
         &         p4(:,:,:,:), &
         &         p5(:,:,:,:,:), &
         &         p6(:,:,:,:,:,:), &
         &         p7(:,:,:,:,:,:,:)

    SELECT CASE (SIZE(coord))

    CASE (1)
      CALL dist_mult_array_cache_ptr_l_1d(dm_array, sub_array_idx, &
           cache_idx, p1)
      v = p1(coord(1))
    CASE (2)
      CALL dist_mult_array_cache_ptr_l_2d(dm_array, sub_array_idx, &
           cache_idx, p2)
      v = p2(coord(1),coord(2))
    CASE (3)
      CALL dist_mult_array_cache_ptr_l_3d(dm_array, sub_array_idx, &
           cache_idx, p3)
      v = p3(coord(1),coord(2),coord(3))
    CASE (4)
      CALL dist_mult_array_cache_ptr_l_4d(dm_array, sub_array_idx, &
           cache_idx, p4)
      v = p4(coord(1),coord(2),coord(3),coord(4))
    CASE (5)
      CALL dist_mult_array_cache_ptr_l_5d(dm_array, sub_array_idx, &
           cache_idx, p5)
      v = p5(coord(1),coord(2),coord(3),coord(4),coord(5))
    CASE (6)
      CALL dist_mult_array_cache_ptr_l_6d(dm_array, sub_array_idx, &
           cache_idx, p6)
      v = p6(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6))
    CASE (7)
      CALL dist_mult_array_cache_ptr_l_7d(dm_array, sub_array_idx, &
           cache_idx, p7)
      v = p7(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6),coord(7))
    CASE default
      CALL abort_ppm("invalid array rank", &
           __FILE__, &
           __LINE__, dm_array%comm)
    END SELECT
  END SUBROUTINE dist_mult_array_cache_val_l


  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_l_1d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    LOGICAL, POINTER :: sub_array_ptr(:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(1)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:1, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_l_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_l_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    LOGICAL, POINTER :: sub_array_ptr(:)

    CALL dist_mult_array_cache_ptr_l_1d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_l_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_l_2d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(2)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:2, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_l_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_l_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:)

    CALL dist_mult_array_cache_ptr_l_2d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_l_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_l_3d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(3)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:3, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_l_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_l_3d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:)

    CALL dist_mult_array_cache_ptr_l_3d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_l_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_l_4d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(4)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:4, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_l_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_l_4d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:)

    CALL dist_mult_array_cache_ptr_l_4d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_l_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_l_5d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(5)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:5, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_l_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_l_5d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_l_5d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_l_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_l_6d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(6)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:6, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_l_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_l_6d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_l_6d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_l_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_l_7d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(7)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:7, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(7, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_l_7d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_l_7d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    LOGICAL, POINTER :: sub_array_ptr(:,:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_l_7d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_l_7d


  SUBROUTINE dist_mult_array_get_deferred_l(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    LOGICAL, INTENT(out) :: v

    INTEGER(mpi_address_kind) :: byte_offset, ofs_factor
    INTEGER :: coord_base(7), a_rank, i
    TYPE(dm_array_cache_entry) :: cache_entry
    INTEGER :: src_comm_rank, comm, dt, ierror
    TYPE(c_ptr) :: baseptr_c
    LOGICAL, POINTER :: baseptr

    comm = dm_array%comm
    src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord)
    ALLOCATE(cache_entry%offset(dm_array%num_sub_arrays))
    CALL compute_cache_addr(cache_entry, dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, src_comm_rank), dm_array%dt_info)
    a_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
    coord_base(1:a_rank) &
         = dm_array%local_chunks(1:a_rank, sub_array, src_comm_rank)%first
    byte_offset = cache_entry%offset(sub_array)
    ofs_factor = dm_array%dt_info(sub_array)%extent
    DO i = 1, a_rank
      byte_offset = byte_offset + (coord(i) - coord_base(i)) * ofs_factor
      ofs_factor = ofs_factor &
           * INT(dm_array%local_chunks(i, sub_array, src_comm_rank)%size, &
           &     mpi_address_kind)
    END DO
    dt = dm_array%sub_arrays_global_desc(sub_array)%element_dt

    IF (src_comm_rank == dm_array%comm_rank) THEN
      baseptr_c = TRANSFER(dm_array%cache(0)%base + byte_offset, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
      CALL mpi_sendrecv(baseptr, 1, dt, 0, 0, v, 1, dt, 0, 0, &
           mpi_comm_self, mpi_status_ignore, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    ELSE
      CALL mpi_get(v, 1, dt, src_comm_rank, byte_offset, 1, dt, dm_array%win, &
           ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    END IF

  END SUBROUTINE dist_mult_array_get_deferred_l

  SUBROUTINE dist_mult_array_get_l(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    LOGICAL, INTENT(out) :: v

    INTEGER :: cache_idx

    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      cache_idx = dist_mult_array_get_cache_idx(dm_array, sub_array, coord, &
           INT(mp_l_extent, mpi_address_kind))
      CALL dist_mult_array_cache_val_l(dm_array, sub_array, &
           cache_idx, coord, v)
    ELSE ! dm_array%sync_mode == sync_mode_active_target
      CALL dist_mult_array_get_deferred_l(dm_array, sub_array, coord, v)
    END IF
  END SUBROUTINE dist_mult_array_get_l


  ! see @ref dist_mult_array_get
  SUBROUTINE dist_mult_array_cache_val_sp(dm_array, sub_array_idx, &
       cache_idx, coord, v)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx, cache_idx
    INTEGER, INTENT(in) :: coord(:)
    REAL(sp), INTENT(out) :: v

    REAL(sp), POINTER :: p1(:), &
         &         p2(:,:), &
         &         p3(:,:,:), &
         &         p4(:,:,:,:), &
         &         p5(:,:,:,:,:), &
         &         p6(:,:,:,:,:,:), &
         &         p7(:,:,:,:,:,:,:)

    SELECT CASE (SIZE(coord))

    CASE (1)
      CALL dist_mult_array_cache_ptr_sp_1d(dm_array, sub_array_idx, &
           cache_idx, p1)
      v = p1(coord(1))
    CASE (2)
      CALL dist_mult_array_cache_ptr_sp_2d(dm_array, sub_array_idx, &
           cache_idx, p2)
      v = p2(coord(1),coord(2))
    CASE (3)
      CALL dist_mult_array_cache_ptr_sp_3d(dm_array, sub_array_idx, &
           cache_idx, p3)
      v = p3(coord(1),coord(2),coord(3))
    CASE (4)
      CALL dist_mult_array_cache_ptr_sp_4d(dm_array, sub_array_idx, &
           cache_idx, p4)
      v = p4(coord(1),coord(2),coord(3),coord(4))
    CASE (5)
      CALL dist_mult_array_cache_ptr_sp_5d(dm_array, sub_array_idx, &
           cache_idx, p5)
      v = p5(coord(1),coord(2),coord(3),coord(4),coord(5))
    CASE (6)
      CALL dist_mult_array_cache_ptr_sp_6d(dm_array, sub_array_idx, &
           cache_idx, p6)
      v = p6(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6))
    CASE (7)
      CALL dist_mult_array_cache_ptr_sp_7d(dm_array, sub_array_idx, &
           cache_idx, p7)
      v = p7(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6),coord(7))
    CASE default
      CALL abort_ppm("invalid array rank", &
           __FILE__, &
           __LINE__, dm_array%comm)
    END SELECT
  END SUBROUTINE dist_mult_array_cache_val_sp


  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_sp_1d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(sp), POINTER :: sub_array_ptr(:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(1)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:1, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_sp_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_sp_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(sp), POINTER :: sub_array_ptr(:)

    CALL dist_mult_array_cache_ptr_sp_1d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_sp_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_sp_2d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(2)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:2, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_sp_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_sp_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:)

    CALL dist_mult_array_cache_ptr_sp_2d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_sp_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_sp_3d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(3)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:3, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_sp_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_sp_3d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:)

    CALL dist_mult_array_cache_ptr_sp_3d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_sp_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_sp_4d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(4)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:4, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_sp_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_sp_4d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:)

    CALL dist_mult_array_cache_ptr_sp_4d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_sp_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_sp_5d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(5)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:5, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_sp_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_sp_5d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_sp_5d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_sp_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_sp_6d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(6)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:6, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_sp_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_sp_6d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_sp_6d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_sp_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_sp_7d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(7)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:7, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(7, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_sp_7d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_sp_7d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(sp), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_sp_7d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_sp_7d


  SUBROUTINE dist_mult_array_get_deferred_sp(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    REAL(sp), INTENT(out) :: v

    INTEGER(mpi_address_kind) :: byte_offset, ofs_factor
    INTEGER :: coord_base(7), a_rank, i
    TYPE(dm_array_cache_entry) :: cache_entry
    INTEGER :: src_comm_rank, comm, dt, ierror
    TYPE(c_ptr) :: baseptr_c
    REAL(sp), POINTER :: baseptr

    comm = dm_array%comm
    src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord)
    ALLOCATE(cache_entry%offset(dm_array%num_sub_arrays))
    CALL compute_cache_addr(cache_entry, dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, src_comm_rank), dm_array%dt_info)
    a_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
    coord_base(1:a_rank) &
         = dm_array%local_chunks(1:a_rank, sub_array, src_comm_rank)%first
    byte_offset = cache_entry%offset(sub_array)
    ofs_factor = dm_array%dt_info(sub_array)%extent
    DO i = 1, a_rank
      byte_offset = byte_offset + (coord(i) - coord_base(i)) * ofs_factor
      ofs_factor = ofs_factor &
           * INT(dm_array%local_chunks(i, sub_array, src_comm_rank)%size, &
           &     mpi_address_kind)
    END DO
    dt = dm_array%sub_arrays_global_desc(sub_array)%element_dt

    IF (src_comm_rank == dm_array%comm_rank) THEN
      baseptr_c = TRANSFER(dm_array%cache(0)%base + byte_offset, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
      CALL mpi_sendrecv(baseptr, 1, dt, 0, 0, v, 1, dt, 0, 0, &
           mpi_comm_self, mpi_status_ignore, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    ELSE
      CALL mpi_get(v, 1, dt, src_comm_rank, byte_offset, 1, dt, dm_array%win, &
           ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    END IF

  END SUBROUTINE dist_mult_array_get_deferred_sp

  SUBROUTINE dist_mult_array_get_sp(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    REAL(sp), INTENT(out) :: v

    INTEGER :: cache_idx

    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      cache_idx = dist_mult_array_get_cache_idx(dm_array, sub_array, coord, &
           INT(mp_sp_extent, mpi_address_kind))
      CALL dist_mult_array_cache_val_sp(dm_array, sub_array, &
           cache_idx, coord, v)
    ELSE ! dm_array%sync_mode == sync_mode_active_target
      CALL dist_mult_array_get_deferred_sp(dm_array, sub_array, coord, v)
    END IF
  END SUBROUTINE dist_mult_array_get_sp


  ! see @ref dist_mult_array_get
  SUBROUTINE dist_mult_array_cache_val_dp(dm_array, sub_array_idx, &
       cache_idx, coord, v)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx, cache_idx
    INTEGER, INTENT(in) :: coord(:)
    REAL(dp), INTENT(out) :: v

    REAL(dp), POINTER :: p1(:), &
         &         p2(:,:), &
         &         p3(:,:,:), &
         &         p4(:,:,:,:), &
         &         p5(:,:,:,:,:), &
         &         p6(:,:,:,:,:,:), &
         &         p7(:,:,:,:,:,:,:)

    SELECT CASE (SIZE(coord))

    CASE (1)
      CALL dist_mult_array_cache_ptr_dp_1d(dm_array, sub_array_idx, &
           cache_idx, p1)
      v = p1(coord(1))
    CASE (2)
      CALL dist_mult_array_cache_ptr_dp_2d(dm_array, sub_array_idx, &
           cache_idx, p2)
      v = p2(coord(1),coord(2))
    CASE (3)
      CALL dist_mult_array_cache_ptr_dp_3d(dm_array, sub_array_idx, &
           cache_idx, p3)
      v = p3(coord(1),coord(2),coord(3))
    CASE (4)
      CALL dist_mult_array_cache_ptr_dp_4d(dm_array, sub_array_idx, &
           cache_idx, p4)
      v = p4(coord(1),coord(2),coord(3),coord(4))
    CASE (5)
      CALL dist_mult_array_cache_ptr_dp_5d(dm_array, sub_array_idx, &
           cache_idx, p5)
      v = p5(coord(1),coord(2),coord(3),coord(4),coord(5))
    CASE (6)
      CALL dist_mult_array_cache_ptr_dp_6d(dm_array, sub_array_idx, &
           cache_idx, p6)
      v = p6(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6))
    CASE (7)
      CALL dist_mult_array_cache_ptr_dp_7d(dm_array, sub_array_idx, &
           cache_idx, p7)
      v = p7(coord(1),coord(2),coord(3),coord(4),coord(5),coord(6),coord(7))
    CASE default
      CALL abort_ppm("invalid array rank", &
           __FILE__, &
           __LINE__, dm_array%comm)
    END SELECT
  END SUBROUTINE dist_mult_array_cache_val_dp


  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_dp_1d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(dp), POINTER :: sub_array_ptr(:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(1)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:1, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_dp_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_dp_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:)

    CALL dist_mult_array_cache_ptr_dp_1d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_dp_1d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_dp_2d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(2)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:2, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_dp_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_dp_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:)

    CALL dist_mult_array_cache_ptr_dp_2d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_dp_2d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_dp_3d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(3)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:3, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_dp_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_dp_3d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:)

    CALL dist_mult_array_cache_ptr_dp_3d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_dp_3d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_dp_4d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(4)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:4, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_dp_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_dp_4d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:)

    CALL dist_mult_array_cache_ptr_dp_4d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_dp_4d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_dp_5d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(5)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:5, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_dp_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_dp_5d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_dp_5d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_dp_5d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_dp_6d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(6)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:6, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_dp_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_dp_6d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_dp_6d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_dp_6d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_cache_ptr_dp_7d(dm_array, sub_array_idx, &
       cache_idx, sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER, INTENT(in) :: cache_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(7)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      __FILE__, &
    !      __LINE__, "invalid data type")
    res_shape = extent_shape(dm_array%local_chunks(1:7, sub_array_idx, &
         dm_array%cache(cache_idx)%rank))
    CALL dist_mult_array_cache_ptr_c(dm_array, sub_array_idx, &
         cache_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr, res_shape)
    sub_array_ptr(extent_start(dm_array%local_chunks(1, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(2, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(3, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(4, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(5, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(6, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):, &
         &        extent_start(dm_array%local_chunks(7, sub_array_idx, &
         &                          dm_array%cache(cache_idx)%rank)):) &
         => sub_array_ptr
  END SUBROUTINE dist_mult_array_cache_ptr_dp_7d

  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE dist_mult_array_local_ptr_dp_7d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:,:,:,:,:,:)

    CALL dist_mult_array_cache_ptr_dp_7d(dm_array, sub_array_idx, 0, &
         sub_array_ptr)
  END SUBROUTINE dist_mult_array_local_ptr_dp_7d


  SUBROUTINE dist_mult_array_get_deferred_dp(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    REAL(dp), INTENT(out) :: v

    INTEGER(mpi_address_kind) :: byte_offset, ofs_factor
    INTEGER :: coord_base(7), a_rank, i
    TYPE(dm_array_cache_entry) :: cache_entry
    INTEGER :: src_comm_rank, comm, dt, ierror
    TYPE(c_ptr) :: baseptr_c
    REAL(dp), POINTER :: baseptr

    comm = dm_array%comm
    src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord)
    ALLOCATE(cache_entry%offset(dm_array%num_sub_arrays))
    CALL compute_cache_addr(cache_entry, dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, src_comm_rank), dm_array%dt_info)
    a_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
    coord_base(1:a_rank) &
         = dm_array%local_chunks(1:a_rank, sub_array, src_comm_rank)%first
    byte_offset = cache_entry%offset(sub_array)
    ofs_factor = dm_array%dt_info(sub_array)%extent
    DO i = 1, a_rank
      byte_offset = byte_offset + (coord(i) - coord_base(i)) * ofs_factor
      ofs_factor = ofs_factor &
           * INT(dm_array%local_chunks(i, sub_array, src_comm_rank)%size, &
           &     mpi_address_kind)
    END DO
    dt = dm_array%sub_arrays_global_desc(sub_array)%element_dt

    IF (src_comm_rank == dm_array%comm_rank) THEN
      baseptr_c = TRANSFER(dm_array%cache(0)%base + byte_offset, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
      CALL mpi_sendrecv(baseptr, 1, dt, 0, 0, v, 1, dt, 0, 0, &
           mpi_comm_self, mpi_status_ignore, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    ELSE
      CALL mpi_get(v, 1, dt, src_comm_rank, byte_offset, 1, dt, dm_array%win, &
           ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
    END IF

  END SUBROUTINE dist_mult_array_get_deferred_dp

#ifdef __BLOCK_GET__
  !NEC: replacement for dist_mult_array_get_deferred_dp, where all data is stored instead of get directly
  SUBROUTINE dist_mult_array_get_dp_blk(dm_array, sub_array, coord, v, proc, ind)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    REAL(dp), INTENT(out) :: v
    INTEGER, INTENT(inout) :: proc 
    INTEGER, INTENT(inout) :: ind 

    INTEGER(mpi_address_kind) :: byte_offset, ofs_factor
    INTEGER :: coord_base(7), a_rank, i
    TYPE(dm_array_cache_entry) :: cache_entry
    INTEGER :: src_comm_rank, comm, dt, ierror
    TYPE(c_ptr) :: baseptr_c
    REAL(dp), POINTER :: baseptr
    INTEGER :: omaxsize
    INTEGER, ALLOCATABLE, DIMENSION(:) :: tmparr

    comm = dm_array%comm
    src_comm_rank = dist_mult_array_coord2rank(dm_array, sub_array, coord)
    ALLOCATE(cache_entry%offset(dm_array%num_sub_arrays))
    CALL compute_cache_addr(cache_entry, dm_array%sub_arrays_global_desc, &
         dm_array%local_chunks(:, :, src_comm_rank), dm_array%dt_info)
    a_rank = dm_array%sub_arrays_global_desc(sub_array)%a_rank
!$NEC novector
    coord_base(1:a_rank) &
         = dm_array%local_chunks(1:a_rank, sub_array, src_comm_rank)%first
    byte_offset = cache_entry%offset(sub_array)
    ofs_factor = dm_array%dt_info(sub_array)%extent
!$NEC novector
    DO i = 1, a_rank
      byte_offset = byte_offset + (coord(i) - coord_base(i)) * ofs_factor
      ofs_factor = ofs_factor &
           * INT(dm_array%local_chunks(i, sub_array, src_comm_rank)%size, &
           &     mpi_address_kind)
    END DO
    dt = dm_array%sub_arrays_global_desc(sub_array)%element_dt

    IF (src_comm_rank == dm_array%comm_rank) THEN
      comm_0_cnt_dp = comm_0_cnt_dp + 1
      baseptr_c = TRANSFER(dm_array%cache(0)%base + byte_offset, baseptr_c)
      CALL C_F_POINTER(baseptr_c, baseptr)
#ifndef __COMM_OPT__
      CALL mpi_sendrecv(baseptr, 1, dt, 0, 0, v, 1, dt, 0, 0, &
           mpi_comm_self, mpi_status_ignore, ierror)
      CALL handle_mpi_error(ierror, comm, __LINE__)
#else
      !NEC: Use assignment instead of MPI-self communication to save overhead
      v = baseptr
#endif
    ELSE
      !NEC: increase access number to this proc
      comm_cnt_dp(src_comm_rank) = comm_cnt_dp(src_comm_rank) + 1
      IF (comm_cnt_dp(src_comm_rank) > comm_lst_dp(src_comm_rank)%maxelem) THEN
        !NEC: reallocation procedure in order to save some memory; tuned, to be efficient
        IF (ALLOCATED(comm_lst_dp(src_comm_rank)%lst)) THEN
          omaxsize = comm_lst_dp(src_comm_rank)%maxelem
          ALLOCATE(tmparr(omaxsize))
          tmparr = comm_lst_dp(src_comm_rank)%lst
          DEALLOCATE(comm_lst_dp(src_comm_rank)%lst)
          comm_lst_dp(src_comm_rank)%maxelem = INT(1.41*omaxsize)
          ALLOCATE(comm_lst_dp(src_comm_rank)%lst(comm_lst_dp(src_comm_rank)%maxelem))
          DO i = 1, omaxsize
            comm_lst_dp(src_comm_rank)%lst(i) = tmparr(i)
          END DO
          DEALLOCATE(tmparr)
        ELSE
          comm_lst_dp(src_comm_rank)%maxelem = lstblksizedp
          ALLOCATE(comm_lst_dp(src_comm_rank)%lst(comm_lst_dp(src_comm_rank)%maxelem))
        END IF
      END IF
     !NEC: store offset and return proc and the access number to this proc
      comm_lst_dp(src_comm_rank)%lst(comm_cnt_dp(src_comm_rank)) = byte_offset
      proc = src_comm_rank
      ind = comm_cnt_dp(src_comm_rank)
    END IF

  END SUBROUTINE dist_mult_array_get_dp_blk

  !NEC: initializes data structure for blocked communication
  SUBROUTINE dist_mult_init_blk_comm_dp(dim1, dim2, dim3, dm_array)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: dim1, dim2, dim3

    !comm_max_dp = dim1 * dim2 * dim3 
    ALLOCATE(comm_cnt_dp(0:dm_array%comm_size-1), comm_lst_dp(0:dm_array%comm_size-1), &
             recv_array_dp(0:dm_array%comm_size-1), index_dbl_array(0:dm_array%comm_size-1))
    comm_lst_dp(:)%maxelem=0
    comm_0_cnt_dp = 0
    comm_cnt_dp = 0
    ALLOCATE(proc_lst_dp(dim1, dim2, dim3), ind_lst_dp(dim1, dim2,dim3))
    proc_lst_dp = -1
  END SUBROUTINE dist_mult_init_blk_comm_dp

  !NEC: executes blocked communcation
  SUBROUTINE dist_mult_do_blk_comm_dp(dm_array)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER :: i, j, ierr
    INTEGER, SAVE :: tsize = 0
   !NEC: Only allocate "ones" once to save time
    ALLOCATE(ones(MAXVAL(comm_cnt_dp)))
    ones = 1
   !NEC: Store MPI type size assuming it is not changing on a process during a run
    IF (tsize == 0) CALL MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, tsize, ierr)
    DO i = 0, dm_array%comm_size-1
      IF (ALLOCATED(comm_lst_dp(i)%lst)) THEN
        !NEC_FU: Make byte offset list into type offset list to be portable
        DO j = 1,comm_cnt_dp(i)
          comm_lst_dp(i)%lst(j) = comm_lst_dp(i)%lst(j) / tsize
        END DO
        call MPI_TYPE_INDEXED(comm_cnt_dp(i), ones, comm_lst_dp(i)%lst, MPI_DOUBLE_PRECISION, index_dbl_array(i), ierr)
        call MPI_TYPE_COMMIT(index_dbl_array(i), ierr)
        ALLOCATE(recv_array_dp(i)%lst(comm_cnt_dp(i)))
        call MPI_GET(recv_array_dp(i)%lst, comm_cnt_dp(i), MPI_DOUBLE_PRECISION, &
                     i, 0_mpi_address_kind, 1, index_dbl_array(i), &
                     dm_array%win, ierr) 
      END IF
    END DO
    DEALLOCATE(ones)
  END SUBROUTINE dist_mult_do_blk_comm_dp

  !NEC: frees data structure for blocked communication
  SUBROUTINE dist_mult_end_blk_comm_dp(dm_array)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER :: i, ierr

    DO i = 0, dm_array%comm_size-1
      IF (ALLOCATED(comm_lst_dp(i)%lst)) THEN
        DEALLOCATE(comm_lst_dp(i)%lst)
        DEALLOCATE(recv_array_dp(i)%lst)
        call MPI_TYPE_FREE(index_dbl_array(i), ierr)
      END IF
    END DO
    !NEC: Make the next selected array size the largest needed, to save
    lstblksizedp = MAX(MAXVAL(comm_cnt_dp),100)
    DEALLOCATE(comm_cnt_dp, comm_lst_dp, recv_array_dp,index_dbl_array)
    DEALLOCATE(proc_lst_dp, ind_lst_dp)
  END SUBROUTINE dist_mult_end_blk_comm_dp
#endif

  SUBROUTINE dist_mult_array_get_dp(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    REAL(dp), INTENT(out) :: v

    INTEGER :: cache_idx

    IF (dm_array%sync_mode == sync_mode_passive_target) THEN
      cache_idx = dist_mult_array_get_cache_idx(dm_array, sub_array, coord, &
           INT(mp_dp_extent, mpi_address_kind))
      CALL dist_mult_array_cache_val_dp(dm_array, sub_array, &
           cache_idx, coord, v)
    ELSE ! dm_array%sync_mode == sync_mode_active_target
      CALL dist_mult_array_get_deferred_dp(dm_array, sub_array, coord, v)
    END IF
  END SUBROUTINE dist_mult_array_get_dp


  SUBROUTINE handle_mpi_error(ierror, comm, line)
    INTEGER, INTENT(in) :: ierror, comm, line
    CHARACTER(256), SAVE :: source = &
__FILE__

    INTEGER :: msg_len, ierror_
    CHARACTER(len=mpi_max_error_string) :: msg

    IF (ierror /= MPI_SUCCESS) THEN
      CALL mpi_error_string(ierror, msg, msg_len, ierror_)
      CALL abort_ppm(msg(1:msg_len), TRIM(source), line, comm)
    END IF
  END SUBROUTINE handle_mpi_error

  SUBROUTINE assertion(cond, source, line, msg)
    LOGICAL, INTENT(in) :: cond
    CHARACTER(*), INTENT(in) :: source, msg
    INTEGER, INTENT(in) :: line
    IF (.NOT. cond) THEN
      CALL abort_ppm("assertion "//msg//" failed", source, line, mpi_comm_world)
    END IF
  END SUBROUTINE assertion

  SUBROUTINE abort_ppm(msg, source, line, comm)
    USE mo_mpi, ONLY: abort_mpi
    CHARACTER(*), INTENT(in) :: source, msg
    INTEGER, INTENT(in) :: line, comm
    WRITE(0, '(4a,i0)') msg, " problem detected at ", source, ":", line
    CALL abort_mpi
  END SUBROUTINE abort_ppm

END MODULE ppm_distributed_array
#else
MODULE ppm_distributed_array
  USE mo_kind, ONLY: i4, i8, sp, dp
  USE mo_exception, ONLY: finish, message_text
  USE ppm_extents, ONLY: extent, is_contained_in
  USE iso_c_binding, ONLY: c_ptr

  INTEGER, PARAMETER :: max_rank = 2

  INTEGER, PARAMETER :: ppm_real_dp = 1
  INTEGER, PARAMETER :: ppm_real_sp = 2
  INTEGER, PARAMETER :: ppm_int = 3
  INTEGER, PARAMETER :: ppm_int_i8 = 4
  INTEGER, PARAMETER :: ppm_bool = 5

  INTEGER, PARAMETER :: not_exposed = 0
  INTEGER, PARAMETER :: exposed = 1

  INTEGER, PARAMETER, PUBLIC :: &
       !> in this mode calls to dist_mult_array_get will immediately
       !! retrieve the requested value
       sync_mode_passive_target = 0, &
       !> in this mode calls to dist_mult_array_get will result in
       !! the passed variable to become defined only after the next call
       !! to dist_mult_array_unexpose
       sync_mode_active_target = 1

  TYPE t_data_ptr
    INTEGER(i4), POINTER :: i4_1d(:)
    INTEGER(i4), POINTER :: i4_2d(:,:)
    INTEGER(i8), POINTER :: i8_1d(:)
    INTEGER(i8), POINTER :: i8_2d(:,:)
    REAL(sp), POINTER :: sp_1d(:)
    REAL(sp), POINTER :: sp_2d(:,:)
    REAL(dp), POINTER :: dp_1d(:)
    REAL(dp), POINTER :: dp_2d(:,:)
  END TYPE t_data_ptr

  TYPE global_array_desc
    INTEGER :: a_rank
    TYPE(extent) :: rect(max_rank)
    INTEGER :: element_dt
  END TYPE global_array_desc

  TYPE dist_mult_array
    PRIVATE
    !> number of arrays that are distributed
    INTEGER :: num_sub_arrays
    !> Per distributed array information on global array shape and contents.\n
    !! The size of this array is 1:num_sub_arrays.
    TYPE(global_array_desc), ALLOCATABLE :: sub_arrays_global_desc(:)
    !> data pointer
    TYPE(t_data_ptr), ALLOCATABLE :: base(:)
    !> exposure status
    INTEGER :: exposure_status
  END TYPE dist_mult_array

  PUBLIC ppm_real_dp, ppm_real_sp, ppm_int, ppm_int_i8, ppm_bool

  PUBLIC :: dist_mult_array, global_array_desc
  PUBLIC :: dist_mult_array_new, dist_mult_array_delete
  PUBLIC :: dist_mult_array_local_ptr, dist_mult_array_get
  PUBLIC :: dist_mult_array_expose, dist_mult_array_unexpose
  PUBLIC :: dist_mult_array_rma_sync

  INTERFACE dist_mult_array_local_ptr
    ! MODULE PROCEDURE dist_mult_array_local_ptr_c
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_1d
    MODULE PROCEDURE dist_mult_array_local_ptr_i4_2d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i4_3d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i4_4d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i4_5d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i4_6d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i4_7d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i8_1d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i8_2d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i8_3d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i8_4d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i8_5d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i8_6d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_i8_7d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_l_1d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_l_2d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_l_3d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_l_4d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_l_5d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_l_6d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_l_7d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_sp_1d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_sp_2d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_sp_3d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_sp_4d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_sp_5d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_sp_6d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_sp_7d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_1d
    MODULE PROCEDURE dist_mult_array_local_ptr_dp_2d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_dp_3d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_dp_4d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_dp_5d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_dp_6d
    ! MODULE PROCEDURE dist_mult_array_local_ptr_dp_7d
  END INTERFACE dist_mult_array_local_ptr

  INTERFACE dist_mult_array_get
    MODULE PROCEDURE dist_mult_array_get_i4
    MODULE PROCEDURE dist_mult_array_get_i8
    ! MODULE PROCEDURE dist_mult_array_get_l
    ! MODULE PROCEDURE dist_mult_array_get_sp
    MODULE PROCEDURE dist_mult_array_get_dp
  END INTERFACE dist_mult_array_get

CONTAINS

  FUNCTION dist_mult_array_new(sub_arrays, local_chunk, comm, cache_size, &
       sync_mode) &
       RESULT(dm_array)
    TYPE(global_array_desc), INTENT(in) :: sub_arrays(:)
    !> shape = (/ max_rank or more, num_sub_arrays /)
    TYPE(extent), INTENT(in) :: local_chunk(:, :)
    INTEGER, INTENT(in) :: comm
    INTEGER, OPTIONAL, INTENT(in) :: cache_size
    INTEGER, OPTIONAL, INTENT(in) :: sync_mode
    TYPE(dist_mult_array) :: dm_array

    INTEGER :: num_sub_arrays, i

    num_sub_arrays = SIZE(sub_arrays)
    dm_array%exposure_status = not_exposed
    dm_array%num_sub_arrays = num_sub_arrays
    ALLOCATE(dm_array%sub_arrays_global_desc(num_sub_arrays), &
      &      dm_array%base(num_sub_arrays))
    dm_array%sub_arrays_global_desc(:) = sub_arrays

    DO i = 1, num_sub_arrays
      IF ((sub_arrays(i)%a_rank == 1) .AND. &
        & (sub_arrays(i)%element_dt == ppm_int)) THEN
        ALLOCATE(dm_array%base(i)%i4_1d( &
          sub_arrays(i)%rect(1)%first:&
          sub_arrays(i)%rect(1)%first + sub_arrays(i)%rect(1)%size - 1))
      ELSE IF (sub_arrays(i)%a_rank == 2 .AND. &
        &      sub_arrays(i)%element_dt == ppm_int) THEN
        ALLOCATE(dm_array%base(i)%i4_2d( &
          sub_arrays(i)%rect(1)%first:sub_arrays(i)%rect(1)%first + &
                                      sub_arrays(i)%rect(1)%size - 1, &
          sub_arrays(i)%rect(2)%first:sub_arrays(i)%rect(2)%first + &
                                      sub_arrays(i)%rect(2)%size - 1))
      ELSE IF (sub_arrays(i)%a_rank == 1 .AND. &
        &      sub_arrays(i)%element_dt == ppm_int_i8) THEN
        ALLOCATE(dm_array%base(i)%i8_1d( &
          sub_arrays(i)%rect(1)%first:&
          sub_arrays(i)%rect(1)%first + sub_arrays(i)%rect(1)%size - 1))
      ELSE IF (sub_arrays(i)%a_rank == 2 .AND. &
        &      sub_arrays(i)%element_dt == ppm_int_i8) THEN
        ALLOCATE(dm_array%base(i)%i8_2d( &
          sub_arrays(i)%rect(1)%first:sub_arrays(i)%rect(1)%first + &
                                      sub_arrays(i)%rect(1)%size - 1, &
          sub_arrays(i)%rect(2)%first:sub_arrays(i)%rect(2)%first + &
                                      sub_arrays(i)%rect(2)%size - 1))
      ELSE IF (sub_arrays(i)%a_rank == 1 .AND. &
        &      sub_arrays(i)%element_dt == ppm_real_sp) THEN
        ALLOCATE(dm_array%base(i)%sp_1d( &
          sub_arrays(i)%rect(1)%first:&
          sub_arrays(i)%rect(1)%first + sub_arrays(i)%rect(1)%size - 1))
      ELSE IF (sub_arrays(i)%a_rank == 2 .AND. &
        &      sub_arrays(i)%element_dt == ppm_real_sp) THEN
        ALLOCATE(dm_array%base(i)%sp_2d( &
          sub_arrays(i)%rect(1)%first:sub_arrays(i)%rect(1)%first + &
                                      sub_arrays(i)%rect(1)%size - 1, &
          sub_arrays(i)%rect(2)%first:sub_arrays(i)%rect(2)%first + &
                                      sub_arrays(i)%rect(2)%size - 1))
      ELSE IF (sub_arrays(i)%a_rank == 1 .AND. &
        &      sub_arrays(i)%element_dt == ppm_real_dp) THEN
        ALLOCATE(dm_array%base(i)%dp_1d( &
          sub_arrays(i)%rect(1)%first:&
          sub_arrays(i)%rect(1)%first + sub_arrays(i)%rect(1)%size - 1))
      ELSE IF (sub_arrays(i)%a_rank == 2 .AND. &
        &      sub_arrays(i)%element_dt == ppm_real_dp) THEN
        ALLOCATE(dm_array%base(i)%dp_2d( &
          sub_arrays(i)%rect(1)%first:sub_arrays(i)%rect(1)%first + &
                                      sub_arrays(i)%rect(1)%size - 1, &
          sub_arrays(i)%rect(2)%first:sub_arrays(i)%rect(2)%first + &
                                      sub_arrays(i)%rect(2)%size - 1))
      END IF
    END DO
  END FUNCTION

  SUBROUTINE dist_mult_array_delete(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER :: i

    DO i = 1, dm_array%num_sub_arrays
      IF (dm_array%sub_arrays_global_desc(i)%a_rank == 1 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_int) THEN
        DEALLOCATE(dm_array%base(i)%i4_1d)
      ELSE IF (dm_array%sub_arrays_global_desc(i)%a_rank == 2 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_int) THEN
        DEALLOCATE(dm_array%base(i)%i4_2d)
      ELSE IF (dm_array%sub_arrays_global_desc(i)%a_rank == 1 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_int_i8) THEN
        DEALLOCATE(dm_array%base(i)%i8_1d)
      ELSE IF (dm_array%sub_arrays_global_desc(i)%a_rank == 2 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_int_i8) THEN
        DEALLOCATE(dm_array%base(i)%i8_2d)
      ELSE IF (dm_array%sub_arrays_global_desc(i)%a_rank == 1 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_real_sp) THEN
        DEALLOCATE(dm_array%base(i)%sp_1d)
      ELSE IF (dm_array%sub_arrays_global_desc(i)%a_rank == 2 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_real_sp) THEN
        DEALLOCATE(dm_array%base(i)%sp_2d)
      ELSE IF (dm_array%sub_arrays_global_desc(i)%a_rank == 1 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_real_dp) THEN
        DEALLOCATE(dm_array%base(i)%dp_1d)
      ELSE IF (dm_array%sub_arrays_global_desc(i)%a_rank == 2 .AND. &
        & dm_array%sub_arrays_global_desc(i)%element_dt == ppm_real_dp) THEN
        DEALLOCATE(dm_array%base(i)%dp_2d)
      END IF
    END DO

    DEALLOCATE(dm_array%sub_arrays_global_desc, dm_array%base)

  END SUBROUTINE dist_mult_array_delete

  SUBROUTINE dist_mult_array_get_i4(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i4), INTENT(out) :: v

    INTEGER :: ref_rank

    CALL assertion(dm_array%exposure_status == exposed, &
         __FILE__, &
         __LINE__, "wrong exposure status")
    CALL assertion(sub_array >= 1 &
         .AND. sub_array <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    ref_rank = SIZE(coord)
    CALL assertion(ref_rank &
         == dm_array%sub_arrays_global_desc(sub_array)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion(is_contained_in(coord, &
         dm_array%sub_arrays_global_desc(sub_array)%rect(1:ref_rank)), &
         __FILE__, &
         __LINE__, "invalid coordinate")

    SELECT CASE (ref_rank)
    CASE(1)
      v = dm_array%base(sub_array)%i4_1d(coord(1))
    CASE(2)
      v = dm_array%base(sub_array)%i4_2d(coord(1), coord(2))
    CASE default
      WRITE(message_text,*) "invalid array rank ", &
        __FILE__, &
        ":", __LINE__
      CALL finish("assertion", message_text)
    END SELECT
  END SUBROUTINE dist_mult_array_get_i4

  SUBROUTINE dist_mult_array_get_i8(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    INTEGER(i8), INTENT(out) :: v

    INTEGER :: ref_rank

    CALL assertion(dm_array%exposure_status == exposed, &
         __FILE__, &
         __LINE__, "wrong exposure status")
    CALL assertion(sub_array >= 1 &
         .AND. sub_array <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    ref_rank = SIZE(coord)
    CALL assertion(ref_rank &
         == dm_array%sub_arrays_global_desc(sub_array)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion(is_contained_in(coord, &
         dm_array%sub_arrays_global_desc(sub_array)%rect(1:ref_rank)), &
         __FILE__, &
         __LINE__, "invalid coordinate")

    SELECT CASE (ref_rank)
    CASE(1)
      v = dm_array%base(sub_array)%i8_1d(coord(1))
    CASE(2)
      v = dm_array%base(sub_array)%i8_2d(coord(1), coord(2))
    CASE default
      WRITE(message_text,*) "invalid array rank ", &
        __FILE__, &
        ":", __LINE__
      CALL finish("assertion", message_text)
    END SELECT
  END SUBROUTINE dist_mult_array_get_i8

  SUBROUTINE dist_mult_array_get_dp(dm_array, sub_array, coord, v)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: sub_array
    INTEGER, INTENT(in) :: coord(:)
    REAL(dp), INTENT(out) :: v

    INTEGER :: ref_rank

    CALL assertion(dm_array%exposure_status == exposed, &
         __FILE__, &
         __LINE__, "wrong exposure status")
    CALL assertion(sub_array >= 1 &
         .AND. sub_array <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    ref_rank = SIZE(coord)
    CALL assertion(ref_rank &
         == dm_array%sub_arrays_global_desc(sub_array)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion(is_contained_in(coord, &
         dm_array%sub_arrays_global_desc(sub_array)%rect(1:ref_rank)), &
         __FILE__, &
         __LINE__, "invalid coordinate")

    SELECT CASE (ref_rank)
    CASE(1)
      v = dm_array%base(sub_array)%dp_1d(coord(1))
    CASE(2)
      v = dm_array%base(sub_array)%dp_2d(coord(1), coord(2))
    CASE default
      WRITE(message_text,*) "invalid array rank ", &
        __FILE__, &
        ":", __LINE__
      CALL finish("assertion", message_text)
    END SELECT
  END SUBROUTINE dist_mult_array_get_dp

  SUBROUTINE dist_mult_array_local_ptr_i4_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:)

    CALL assertion(sub_array_idx >= 1 &
         .AND. sub_array_idx <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    CALL assertion( &
         1 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion( ppm_int &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         __FILE__, &
         __LINE__, "type mismatch in array reference")

    sub_array_ptr => dm_array%base(sub_array_idx)%i4_1d(:)

  END SUBROUTINE dist_mult_array_local_ptr_i4_1d

  SUBROUTINE dist_mult_array_local_ptr_i4_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:)

    CALL assertion(sub_array_idx >= 1 &
         .AND. sub_array_idx <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    CALL assertion( &
         2 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion( ppm_int &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         __FILE__, &
         __LINE__, "type mismatch in array reference")

    sub_array_ptr => dm_array%base(sub_array_idx)%i4_2d(:,:)

  END SUBROUTINE dist_mult_array_local_ptr_i4_2d

  SUBROUTINE dist_mult_array_local_ptr_dp_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:)

    CALL assertion(sub_array_idx >= 1 &
         .AND. sub_array_idx <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    CALL assertion( &
         1 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion( ppm_real_dp &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         __FILE__, &
         __LINE__, "type mismatch in array reference")

    sub_array_ptr => dm_array%base(sub_array_idx)%dp_1d(:)

  END SUBROUTINE dist_mult_array_local_ptr_dp_1d

  SUBROUTINE dist_mult_array_local_ptr_dp_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:)

    CALL assertion(sub_array_idx >= 1 &
         .AND. sub_array_idx <= SIZE(dm_array%sub_arrays_global_desc), &
         __FILE__, &
         __LINE__, "invalid subarray index")
    CALL assertion( &
         2 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         __FILE__, &
         __LINE__, "rank mismatch in array reference")
    CALL assertion( ppm_real_dp &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         __FILE__, &
         __LINE__, "type mismatch in array reference")

    sub_array_ptr => dm_array%base(sub_array_idx)%dp_2d

  END SUBROUTINE dist_mult_array_local_ptr_dp_2d

  SUBROUTINE dist_mult_array_expose(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    dm_array%exposure_status = exposed
  END SUBROUTINE dist_mult_array_expose

  SUBROUTINE dist_mult_array_unexpose(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    dm_array%exposure_status = not_exposed
  END SUBROUTINE dist_mult_array_unexpose

  SUBROUTINE dist_mult_array_rma_sync(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
  END SUBROUTINE dist_mult_array_rma_sync

  SUBROUTINE assertion(cond, source, line, msg)
    LOGICAL, INTENT(in) :: cond
    CHARACTER(*), INTENT(in) :: source, msg
    INTEGER, INTENT(in) :: line
    IF (.NOT. cond) THEN
      WRITE(message_text,'(5a,i0)') "assertion ", msg, " failed ", source, &
           ":", line
      CALL finish("assertion", message_text)
    END IF
  END SUBROUTINE assertion
END MODULE
#endif
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! mode: f90
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
