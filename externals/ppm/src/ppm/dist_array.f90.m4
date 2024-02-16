changequote(`{',`}')dnl
include({forloop2.m4})dnl
!> @file dist_array.f90
!! @brief distributed array implementation,
!! optimized for read-only access in parallel phase
!!
!!
!! @copyright Copyright  (C)  2016  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
define({interface_gen},
{forloop({dim},{1},{7},{{    MODULE PROCEDURE $2_$1_}dim{d}
})})

MODULE ppm_distributed_array_rename
  IMPLICIT NONE
  PUBLIC
  EXTERNAL :: ppm_dist_mult_array_get_f2c
END MODULE ppm_distributed_array_rename

MODULE ppm_distributed_array
  !USE ppm_base, ONLY: assertion
  USE ppm_std_type_kinds, ONLY: i4, i8, sp, dp
  USE ppm_extents, ONLY: extent, extent_shape, &
       extent_start
  USE iso_c_binding, ONLY: c_ptr, c_f_pointer
  USE ppm_distributed_array_rename, ONLY: &
       dist_mult_array_get => ppm_dist_mult_array_get_f2c
define({ptr_bnds_remap},{{
#ifdef HAVE_POINTER_REMAP
    $1($3(1):}forloop({pdim},{2},dim,{, $3(pdim):}){) &
         => $2
#else
    CALL ptr_bnds_remap($1, $2, $3)
#endif
}})dnl
#ifndef HAVE_POINTER_REMAP
  {USE ppm_ptr_bnds_remap, ONLY: ptr_bnds_remap}
#endif
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  INCLUDE 'ftype_size.inc'
  INTEGER, PARAMETER :: max_rank = 7
  !> describes one array-like global data structure to be distributed
  TYPE global_array_desc
    SEQUENCE
    !> Fortran array rank of array
    INTEGER :: a_rank
    !> rect(1:a_rank) describes bounds of global array
    TYPE(extent) :: rect(max_rank)
    !> MPI datatype that describes the elements of an array
    INTEGER :: element_dt
  END TYPE global_array_desc

  INTEGER, PARAMETER, PUBLIC :: &
       !> in this mode calls to dist_mult_array_get will immediately
       !! retrieve the requested value
       sync_mode_passive_target = 0, &
       !> in this mode calls to dist_mult_array_get will result in
       !! the passed variable to become defined only after the next call
       !! to dist_mult_array_unexpose
       sync_mode_active_target = 1, &
       !> in this mode, remote access is disabled at first but can later be
       !! activated
       sync_mode_local_only = 2

  INTEGER, PARAMETER, PUBLIC :: &
       transfer_mode_struct = 0, &
       transfer_mode_bytes = 1

  !> dist_mult_array describes a global array where each rank only
  !! holds a part of the whole and this object can be used to access values
  !! stored on other ranks.
  TYPE dist_mult_array
    PRIVATE
    SEQUENCE
    INTEGER(mpi_address_kind) :: handle
  END TYPE dist_mult_array

  PUBLIC :: dist_mult_array, global_array_desc
  PUBLIC :: dist_mult_array_new, dist_mult_array_delete
  PUBLIC :: dist_mult_array_copy
  PUBLIC :: dist_mult_array_local_ptr, dist_mult_array_get
  PUBLIC :: dist_mult_array_expose, dist_mult_array_unexpose
  PUBLIC :: dist_mult_array_rma_sync
  PUBLIC :: dist_mult_array_comm
  PUBLIC :: dist_mult_array_set_transfer_mode, dist_mult_array_get_transfer_mode
  PUBLIC :: dist_mult_array_get_sync_mode
  PUBLIC :: dist_mult_array_set_sync_mode

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
interface_gen({i4},{dist_mult_array_local_ptr})dnl
interface_gen({i8},{dist_mult_array_local_ptr})dnl
interface_gen({l},{dist_mult_array_local_ptr})dnl
interface_gen({sp},{dist_mult_array_local_ptr})dnl
interface_gen({dp},{dist_mult_array_local_ptr})dnl
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
  !!

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

    INTEGER :: num_sub_arrays, cache_size_, max_a_rank, sync_mode_
    TYPE(extent), ALLOCATABLE :: local_chunk_(:,:)
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_new_f2c(new_dma_handle, num_sub_arrays, &
           sub_arrays, local_chunk, comm, cache_size, sync_mode)
        IMPORT :: mpi_address_kind, global_array_desc, extent
        INTEGER(mpi_address_kind), INTENT(out) :: new_dma_handle
        INTEGER, INTENT(in) :: num_sub_arrays, comm, cache_size, sync_mode
        TYPE(global_array_desc), INTENT(in) :: sub_arrays(*)
        TYPE(extent), INTENT(in) :: local_chunk(*)
      END SUBROUTINE ppm_dist_mult_array_new_f2c
    END INTERFACE

    num_sub_arrays = SIZE(sub_arrays)
    IF (PRESENT(cache_size)) THEN ; cache_size_ = cache_size ; ELSE
      cache_size_ = 0
    END IF
    IF (PRESENT(sync_mode)) THEN ; sync_mode_ = sync_mode ; ELSE
      sync_mode_ =  sync_mode_passive_target
    END IF
    max_a_rank = MAXVAL(sub_arrays%a_rank)
    ALLOCATE(local_chunk_(max_a_rank, num_sub_arrays))
    local_chunk_ = local_chunk(1:max_a_rank, 1:num_sub_arrays)
    CALL ppm_dist_mult_array_new_f2c(dm_array%handle, num_sub_arrays, &
         sub_arrays, local_chunk_, comm, cache_size_, sync_mode_)
  END FUNCTION dist_mult_array_new

  !> @brief copy dist_mult_array data type
  !! @remark since an MPI window is created, this routine is collective for
  !! all processes which participated in the corresponding call to
  !! @ref dist_mult_array_new or @ref dist_mult_array_copy
  !! @param[inout] dm_array distributed multi-array to copy
  !! @return copy of distributed multi-array \a dm_array
  FUNCTION dist_mult_array_copy(dm_array) RESULT(dm_array_copy)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    TYPE(dist_mult_array) :: dm_array_copy
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_copy_f2c(handle, copy_handle)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: handle
        INTEGER(mpi_address_kind), INTENT(out) :: copy_handle
      END SUBROUTINE ppm_dist_mult_array_copy_f2c
    END INTERFACE
    CALL ppm_dist_mult_array_copy_f2c(dm_array%handle, dm_array_copy%handle)
  END FUNCTION dist_mult_array_copy

  !> @brief destruct dist_mult_array data type
  !! @remark since an MPI window is freed, this routine is collective for
  !! all processes which participated in the corresponding call to
  !! @ref dist_mult_array_new or @ref dist_mult_array_copy
  !! @param[inout] dm_array distributed multi-array to delete
  SUBROUTINE dist_mult_array_delete(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_delete_f2c(handle)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: handle
      END SUBROUTINE ppm_dist_mult_array_delete_f2c
    END INTERFACE
    CALL ppm_dist_mult_array_delete_f2c(dm_array%handle)
  END SUBROUTINE dist_mult_array_delete

  !> retrieve TYPE(c_ptr) to local chunk
  SUBROUTINE dist_mult_array_local_ptr_c(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    TYPE(c_ptr), INTENT(out) :: sub_array_ptr
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_local_ptr_f2c(handle, sub_array_idx, cptr)
        IMPORT :: c_ptr, mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: handle
        INTEGER, INTENT(in) :: sub_array_idx
        TYPE(c_ptr), INTENT(out) :: cptr
      END SUBROUTINE ppm_dist_mult_array_local_ptr_f2c
    END INTERFACE

    CALL ppm_dist_mult_array_local_ptr_f2c(dm_array%handle, sub_array_idx, &
         sub_array_ptr)
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

    INTERFACE
      SUBROUTINE ppm_dist_mult_array_expose_f2c(handle)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: handle
      END SUBROUTINE ppm_dist_mult_array_expose_f2c
    END INTERFACE
    CALL ppm_dist_mult_array_expose_f2c(dm_array%handle)
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

    INTERFACE
      SUBROUTINE ppm_dist_mult_array_unexpose_f2c(handle)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: handle
      END SUBROUTINE ppm_dist_mult_array_unexpose_f2c
    END INTERFACE
    CALL ppm_dist_mult_array_unexpose_f2c(dm_array%handle)
  END SUBROUTINE dist_mult_array_unexpose

  !> @brief synchronize RMA updates only, ignore local updates
  !!
  !! This call is collective for all ranks in the communicator used
  !! to create @a dm_array.
  !! Also @a dm_array must be in exposed state.
  SUBROUTINE dist_mult_array_rma_sync(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_rma_sync_f2c(handle)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: handle
      END SUBROUTINE ppm_dist_mult_array_rma_sync_f2c
    END INTERFACE
    CALL ppm_dist_mult_array_rma_sync_f2c(dm_array%handle)
  END SUBROUTINE dist_mult_array_rma_sync

  !> @brief query sync protocol
  !!
  FUNCTION dist_mult_array_get_sync_mode(dm_array) RESULT(sync_mode)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER :: sync_mode
    INTERFACE
      FUNCTION PPM_dist_mult_array_get_sync_mode_f2c(handle) RESULT(sync_mode)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: handle
        INTEGER :: sync_mode
      END FUNCTION PPM_dist_mult_array_get_sync_mode_f2c
    END INTERFACE
    sync_mode = ppm_dist_mult_array_get_sync_mode_f2c(dm_array%handle)
  END FUNCTION dist_mult_array_get_sync_mode

  !> @brief change sync protocol
  !!
  !! This call is collective for all ranks in the communicator used
  !! to create @a dm_array
  !!
  !! note: this ends an exposure epoch on the data shared via RMA
  SUBROUTINE dist_mult_array_set_sync_mode(dm_array, sync_mode, cache_size)
    !> distributed multiple arrays to change access/sync mode for
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    !> new sync mode to set
    INTEGER, INTENT(in) :: sync_mode
    !> in case @a sync_mode is sync_mode_passive_target, the number of
    !! ranks to cache data for is set to this value (or determined
    !! automatically if 0 or not present)
    INTEGER, OPTIONAL, INTENT(in) :: cache_size
    INTEGER :: cache_size_
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_set_sync_mode_f2c(handle, sync_mode, &
           cache_size)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind) :: handle
        INTEGER, INTENT(in) :: sync_mode, cache_size
      END SUBROUTINE ppm_dist_mult_array_set_sync_mode_f2c
    END INTERFACE
    IF (PRESENT(cache_size)) THEN
      cache_size_ = cache_size
    ELSE
      cache_size_ = 0
    END IF
    CALL ppm_dist_mult_array_set_sync_mode_f2c(dm_array%handle, sync_mode, &
         cache_size_)
  END SUBROUTINE dist_mult_array_set_sync_mode

  SUBROUTINE dist_mult_array_rank_rect(dm_array, sub_array_idx, rank, rect)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx, rank
    TYPE(extent), INTENT(out) :: rect(*)
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_rank_rect_f2c(dma_handle, sub_array_idx, &
           rank, rect)
        IMPORT :: mpi_address_kind, extent
        INTEGER(mpi_address_kind), INTENT(in) :: dma_handle
        INTEGER, INTENT(in) :: sub_array_idx, rank
        TYPE(extent), INTENT(out) :: rect(*)
      END SUBROUTINE PPM_dist_mult_array_rank_rect_f2c
    END INTERFACE
    CALL ppm_dist_mult_array_rank_rect_f2c(dm_array%handle, sub_array_idx, &
         rank, rect)
  END SUBROUTINE dist_mult_array_rank_rect

  FUNCTION dist_mult_array_comm(dm_array) RESULT(comm)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER :: comm
    INTERFACE
      FUNCTION ppm_dist_mult_array_comm_f2c(dma_handle) RESULT(comm_f)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: dma_handle
        INTEGER :: comm_f
      END FUNCTION ppm_dist_mult_array_comm_f2c
    END INTERFACE
    comm = ppm_dist_mult_array_comm_f2c(dm_array%handle)
  END FUNCTION dist_mult_array_comm

  SUBROUTINE dist_mult_array_set_transfer_mode(dm_array, mode)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER, INTENT(in) :: mode
    INTERFACE
      SUBROUTINE ppm_dist_mult_array_set_transfer_mode_f2c(dma_handle, mode)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(inout) :: dma_handle
        INTEGER, INTENT(in) :: mode
      END SUBROUTINE ppm_dist_mult_array_set_transfer_mode_f2c
    END INTERFACE
    CALL ppm_dist_mult_array_set_transfer_mode_f2c(dm_array%handle, mode)
  END SUBROUTINE dist_mult_array_set_transfer_mode

  FUNCTION dist_mult_array_get_transfer_mode(dm_array) RESULT(mode)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array
    INTEGER :: mode
    INTERFACE
      FUNCTION ppm_dist_mult_array_get_transfer_mode_f2c(dma_handle) &
           RESULT(mode_c)
        IMPORT :: mpi_address_kind
        INTEGER(mpi_address_kind), INTENT(in) :: dma_handle
        INTEGER :: mode_c
      END FUNCTION ppm_dist_mult_array_get_transfer_mode_f2c
    END INTERFACE
    mode = ppm_dist_mult_array_get_transfer_mode_f2c(dm_array%handle)
  END FUNCTION dist_mult_array_get_transfer_mode

define({type_gen},{
forloop({dim},{1},{7},{
  ! see @ref dist_mult_array_local_ptr
  SUBROUTINE {dist_mult_array_local_ptr_$1_}dim{d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    $2, POINTER, DIMENSION(:}forloop({pdim},{2},dim,{,:}){) :: sub_array_ptr, &
         sub_array_ptr1
    TYPE(c_ptr) :: a_ptr
    INTEGER :: res_shape(}dim{)
    TYPE(extent) :: res_rect(}dim{)
    INTEGER :: lb(}dim{)
    ! such a check as follows would be nice here, but is impossible with
    ! current mpi standard
    ! CALL assertion(dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt &
    !      == mp_i4, &
    !      filename, &
    !      __LINE__, "invalid data type")
    CALL dist_mult_array_rank_rect(dm_array, sub_array_idx, -1, res_rect)
    res_shape = extent_shape(res_rect)
    CALL dist_mult_array_local_ptr_c(dm_array, sub_array_idx, a_ptr)
    CALL C_F_POINTER(a_ptr, sub_array_ptr1, res_shape)
    lb = extent_start(res_rect)}ptr_bnds_remap(sub_array_ptr, sub_array_ptr1, lb){
  END SUBROUTINE dist_mult_array_local_ptr_$1_}dim{d
}})})
type_gen({i4},{INTEGER(i4)})
type_gen({i8},{INTEGER(i8)})
type_gen({l},{LOGICAL})
type_gen({sp},{REAL(sp)})
type_gen({dp},{REAL(dp)})

END MODULE ppm_distributed_array
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! mode: f90
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
