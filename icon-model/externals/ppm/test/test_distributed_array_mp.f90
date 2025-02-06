!> @file test_distributed_array_mp.f90
!! @brief test of distributed array
!!
!! @copyright Copyright  (C)  2013  Thomas Jahns <jahns@dkrz.de>
!!
!! @version 1.0
!! @author Thomas Jahns <jahns@dkrz.de>
!
! Keywords:
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
#include "fc_feature_defs.inc"
PROGRAM test_dist_array
  USE ppm_base, ONLY: abort_ppm, assertion
  USE ppm_distributed_array, ONLY: dist_mult_array, global_array_desc, &
       dist_mult_array_new, dist_mult_array_delete, &
       dist_mult_array_local_ptr, dist_mult_array_get, &
       dist_mult_array_expose, dist_mult_array_unexpose, &
       dist_mult_array_get_sync_mode, dist_mult_array_set_sync_mode, &
       sync_mode_active_target, sync_mode_passive_target, sync_mode_local_only
  USE ppm_extents, ONLY: extent
  USE ppm_std_type_kinds, ONLY: i4
  USE scales_ppm, ONLY: initialize_scales_ppm, finalize_scales_ppm
#ifdef USE_MPI_MOD
  USE mpi
#endif
  IMPLICIT NONE
#if defined USE_MPI && ! defined USE_MPI_MOD
  INCLUDE 'mpif.h'
#endif
  TYPE(dist_mult_array) :: dm_array
  INTEGER, PARAMETER :: num_arrays = 2, max_rank = 2
  TYPE(global_array_desc) :: glob_spec(num_arrays)
  TYPE(extent) :: local_array(max_rank, num_arrays)
  INTEGER :: ierror, comm_rank, comm_size
  INTEGER(i4) :: random_seed
  INTEGER :: cache_mode
  INTEGER, PARAMETER :: cache_mode_auto = 1, cache_mode_all = 2
  CHARACTER(len=*), PARAMETER :: filename = 'test_distributed_array_mp.f90'
  ! init mpi
  CALL mpi_init(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_init failed", filename, __LINE__)
  ! init scales ppm
  CALL initialize_scales_ppm(random_seed=0, seed_output=random_seed)
  ! find rank and number of tasks
  CALL mpi_comm_rank(mpi_comm_world, comm_rank, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_comm_rank failed", filename, __LINE__)
  CALL mpi_comm_size(mpi_comm_world, comm_size, ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_comm_size failed", filename, __LINE__)
  PRINT '(i0,a,i0)', comm_rank, ': random_seed=', random_seed
  ! initialize distributed array
  glob_spec(1) = global_array_desc(1, extent(1, comm_size), mpi_integer)
  glob_spec(2) = global_array_desc(2, extent(1, 2), &
       mpi_double_precision)
  glob_spec(2)%rect(2) = extent(0, comm_size)
  local_array(1, 1) = extent(comm_rank + 1, 1)
  local_array(1, 2) = extent(1, 2)
  local_array(2, 2) = extent(comm_rank, 1)
  dm_array = dist_mult_array_new(glob_spec, local_array, mpi_comm_world, &
       sync_mode=sync_mode_local_only)
  CALL test_array_data(dm_array)
  CALL dist_mult_array_delete(dm_array)
  DO cache_mode = 1, 2
    IF (cache_mode == cache_mode_auto) THEN
      ! run with default caching at first ...
      dm_array = dist_mult_array_new(glob_spec, local_array, mpi_comm_world)
    ELSE
      ! ... and no cache-evictions second
      dm_array = dist_mult_array_new(glob_spec, local_array, mpi_comm_world, &
           comm_size)
    END IF
    CALL test_array_data(dm_array)
    CALL dist_mult_array_delete(dm_array)
  END DO
  dm_array = dist_mult_array_new(glob_spec, local_array, mpi_comm_world, &
       sync_mode=sync_mode_active_target)
  CALL test_array_data(dm_array)
  CALL dist_mult_array_delete(dm_array)

  ! next test switching from one mode to another
  ! start test with local data only
  dm_array = dist_mult_array_new(glob_spec, local_array, mpi_comm_world, &
       sync_mode=sync_mode_local_only)
  CALL test_array_data(dm_array)
  ! switch local -> passive target
  DO cache_mode = 1, 2
    ! run with little caching at first ...
    ! ... and no cache-evictions second
    CALL dist_mult_array_set_sync_mode(dm_array, sync_mode_passive_target, &
         MERGE(0, comm_size, cache_mode == 1))
    CALL test_array_data(dm_array)
  END DO
  ! switch passive -> active
  CALL dist_mult_array_set_sync_mode(dm_array, sync_mode_active_target)
  CALL test_array_data(dm_array)
  ! switch active -> local
  CALL dist_mult_array_set_sync_mode(dm_array, sync_mode_local_only)
  CALL test_array_data(dm_array)
  ! switch local -> active
  CALL dist_mult_array_set_sync_mode(dm_array, sync_mode_active_target)
  CALL test_array_data(dm_array)
  ! switch active -> passive
  DO cache_mode = 1, 2
    ! run with little caching at first ...
    ! ... and no cache-evictions second
    CALL dist_mult_array_set_sync_mode(dm_array, sync_mode_passive_target, &
         MERGE(0, comm_size, cache_mode == 1))
    CALL test_array_data(dm_array)
  END DO
  ! switch passive -> local
  CALL dist_mult_array_set_sync_mode(dm_array, sync_mode_local_only)
  CALL test_array_data(dm_array)
  CALL dist_mult_array_delete(dm_array)

  CALL finalize_scales_ppm
  CALL mpi_finalize(ierror)
  IF (ierror /= mpi_success) &
       CALL abort_ppm("mpi_finalize failed", filename, __LINE__)
CONTAINS
  SUBROUTINE test_array_data(dm_array)
    TYPE(dist_mult_array), INTENT(inout) :: dm_array

    DOUBLE PRECISION, PARAMETER :: eps = 1.0d-13
    DOUBLE PRECISION, ALLOCATABLE :: values_d(:, :)
    DOUBLE PRECISION, POINTER :: local_d_ref(:, :)
    INTEGER, ALLOCATABLE :: values_i(:)
    INTEGER, POINTER :: local_i_ref(:)
    INTEGER :: i, glob_idx, coord(2)
    INTEGER :: sync_mode

    sync_mode = dist_mult_array_get_sync_mode(dm_array)
    ! fill local portion of distributed array
    CALL dist_mult_array_local_ptr(dm_array, 1, local_i_ref)
    local_i_ref(comm_rank + 1) = comm_rank
    ! expose local part available
    CALL dist_mult_array_expose(dm_array)
    ! get values from other ranks
    ALLOCATE(values_i(comm_size))
    IF (sync_mode == sync_mode_passive_target) THEN
      DO glob_idx = 1, comm_size
        coord(2) = glob_idx
        CALL dist_mult_array_get(dm_array, 1, coord(2:2), values_i(glob_idx))
        CALL assertion(values_i(glob_idx) == glob_idx - 1, filename, __LINE__, &
             "Unexpected value in array.")
      END DO
      CALL dist_mult_array_unexpose(dm_array)
    ELSE IF (sync_mode == sync_mode_active_target) THEN
      DO glob_idx = 1, comm_size
        coord(2) = glob_idx
        CALL dist_mult_array_get(dm_array, 1, coord(2:2), values_i(glob_idx))
      END DO
      CALL dist_mult_array_unexpose(dm_array)
      DO glob_idx = 0, comm_size - 1
        CALL assertion(values_i(glob_idx + 1) == glob_idx, filename, __LINE__, &
             "Unexpected value in array.")
      END DO
    ELSE IF (sync_mode == sync_mode_local_only) THEN
      coord(2) = comm_rank+1
      CALL dist_mult_array_get(dm_array, 1, coord(2:2), values_i(1))
      CALL assertion(values_i(1) == comm_rank, filename, __LINE__, &
           "Unexpected value in array.")
      CALL dist_mult_array_unexpose(dm_array)
    END IF
    ! second phase: also use double array
    CALL dist_mult_array_local_ptr(dm_array, 2, local_d_ref)
    local_d_ref(1, comm_rank) = LOG10(DBLE(comm_rank + 3))
    local_d_ref(2, comm_rank) = EXP(DBLE(comm_rank))
    CALL dist_mult_array_expose(dm_array)
    ALLOCATE(values_d(2, 0 : comm_size - 1))
    IF (sync_mode == sync_mode_passive_target) THEN
      DO glob_idx = 0, comm_size - 1
        coord(2) = glob_idx
        DO i = 1, 2
          coord(1) = i
          CALL dist_mult_array_get(dm_array, 2, coord, values_d(i, glob_idx))
        END DO
        CALL assertion(&
             ABS(values_d(1, glob_idx) - LOG10(DBLE(glob_idx + 3))) < eps &
             .AND. ABS(values_d(2, glob_idx) - EXP(DBLE(glob_idx))) < eps, &
             filename, __LINE__, "Unexpected value in array.")
      END DO
    ELSE IF (sync_mode == sync_mode_active_target) THEN
      DO glob_idx = 0, comm_size - 1
        coord(2) = glob_idx
        DO i = 1, 2
          coord(1) = i
          CALL dist_mult_array_get(dm_array, 2, coord, values_d(i, glob_idx))
        END DO
      END DO
      CALL dist_mult_array_unexpose(dm_array)
      DO glob_idx = 0, comm_size - 1
        CALL assertion(&
             ABS(values_d(1, glob_idx) - LOG10(DBLE(glob_idx + 3))) < eps &
             .AND. ABS(values_d(2, glob_idx) - EXP(DBLE(glob_idx))) < eps, &
             filename, __LINE__, "Unexpected value in array.")
      END DO
    ELSE IF (sync_mode == sync_mode_local_only) THEN
      coord(2) = comm_rank
      DO i = 1, 2
        coord(1) = i
        CALL dist_mult_array_get(dm_array, 2, coord, values_d(i, 0))
      END DO
      CALL assertion(&
           ABS(values_d(1, 0) - LOG10(DBLE(comm_rank + 3))) < eps &
           .AND. ABS(values_d(2, 0) - EXP(DBLE(comm_rank))) < eps, &
           filename, __LINE__, "Unexpected value in array.")
    END IF
    DEALLOCATE(values_d, values_i)
  END SUBROUTINE test_array_data
END PROGRAM test_dist_array
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
