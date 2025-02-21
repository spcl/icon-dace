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
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3489, "wrong exposure status")
    CALL assertion(sub_array >= 1 &
         .AND. sub_array <= SIZE(dm_array%sub_arrays_global_desc), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3493, "invalid subarray index")
    ref_rank = SIZE(coord)
    CALL assertion(ref_rank &
         == dm_array%sub_arrays_global_desc(sub_array)%a_rank, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3498, "rank mismatch in array reference")
    CALL assertion(is_contained_in(coord, &
         dm_array%sub_arrays_global_desc(sub_array)%rect(1:ref_rank)), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3502, "invalid coordinate")

    SELECT CASE (ref_rank)
    CASE(1)
      v = dm_array%base(sub_array)%i4_1d(coord(1))
    CASE(2)
      v = dm_array%base(sub_array)%i4_2d(coord(1), coord(2))
    CASE default
      WRITE(message_text,*) "invalid array rank ", &
        "src/parallel_infrastructure/mo_dist_array.f90", &
        ":", 3512
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
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3527, "wrong exposure status")
    CALL assertion(sub_array >= 1 &
         .AND. sub_array <= SIZE(dm_array%sub_arrays_global_desc), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3531, "invalid subarray index")
    ref_rank = SIZE(coord)
    CALL assertion(ref_rank &
         == dm_array%sub_arrays_global_desc(sub_array)%a_rank, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3536, "rank mismatch in array reference")
    CALL assertion(is_contained_in(coord, &
         dm_array%sub_arrays_global_desc(sub_array)%rect(1:ref_rank)), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3540, "invalid coordinate")

    SELECT CASE (ref_rank)
    CASE(1)
      v = dm_array%base(sub_array)%i8_1d(coord(1))
    CASE(2)
      v = dm_array%base(sub_array)%i8_2d(coord(1), coord(2))
    CASE default
      WRITE(message_text,*) "invalid array rank ", &
        "src/parallel_infrastructure/mo_dist_array.f90", &
        ":", 3550
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
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3565, "wrong exposure status")
    CALL assertion(sub_array >= 1 &
         .AND. sub_array <= SIZE(dm_array%sub_arrays_global_desc), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3569, "invalid subarray index")
    ref_rank = SIZE(coord)
    CALL assertion(ref_rank &
         == dm_array%sub_arrays_global_desc(sub_array)%a_rank, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3574, "rank mismatch in array reference")
    CALL assertion(is_contained_in(coord, &
         dm_array%sub_arrays_global_desc(sub_array)%rect(1:ref_rank)), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3578, "invalid coordinate")

    SELECT CASE (ref_rank)
    CASE(1)
      v = dm_array%base(sub_array)%dp_1d(coord(1))
    CASE(2)
      v = dm_array%base(sub_array)%dp_2d(coord(1), coord(2))
    CASE default
      WRITE(message_text,*) "invalid array rank ", &
        "src/parallel_infrastructure/mo_dist_array.f90", &
        ":", 3588
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
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3602, "invalid subarray index")
    CALL assertion( &
         1 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3606, "rank mismatch in array reference")
    CALL assertion( ppm_int &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3610, "type mismatch in array reference")

    sub_array_ptr => dm_array%base(sub_array_idx)%i4_1d(:)

  END SUBROUTINE dist_mult_array_local_ptr_i4_1d

  SUBROUTINE dist_mult_array_local_ptr_i4_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    INTEGER(i4), POINTER :: sub_array_ptr(:,:)

    CALL assertion(sub_array_idx >= 1 &
         .AND. sub_array_idx <= SIZE(dm_array%sub_arrays_global_desc), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3625, "invalid subarray index")
    CALL assertion( &
         2 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3629, "rank mismatch in array reference")
    CALL assertion( ppm_int &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3633, "type mismatch in array reference")

    sub_array_ptr => dm_array%base(sub_array_idx)%i4_2d(:,:)

  END SUBROUTINE dist_mult_array_local_ptr_i4_2d

  SUBROUTINE dist_mult_array_local_ptr_dp_1d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:)

    CALL assertion(sub_array_idx >= 1 &
         .AND. sub_array_idx <= SIZE(dm_array%sub_arrays_global_desc), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3648, "invalid subarray index")
    CALL assertion( &
         1 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3652, "rank mismatch in array reference")
    CALL assertion( ppm_real_dp &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3656, "type mismatch in array reference")

    sub_array_ptr => dm_array%base(sub_array_idx)%dp_1d(:)

  END SUBROUTINE dist_mult_array_local_ptr_dp_1d

  SUBROUTINE dist_mult_array_local_ptr_dp_2d(dm_array, sub_array_idx, &
       sub_array_ptr)
    TYPE(dist_mult_array), INTENT(in) :: dm_array
    INTEGER, INTENT(in) :: sub_array_idx
    REAL(dp), POINTER :: sub_array_ptr(:,:)

    CALL assertion(sub_array_idx >= 1 &
         .AND. sub_array_idx <= SIZE(dm_array%sub_arrays_global_desc), &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3671, "invalid subarray index")
    CALL assertion( &
         2 == dm_array%sub_arrays_global_desc(sub_array_idx)%a_rank, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3675, "rank mismatch in array reference")
    CALL assertion( ppm_real_dp &
         == dm_array%sub_arrays_global_desc(sub_array_idx)%element_dt, &
         "src/parallel_infrastructure/mo_dist_array.f90", &
         3679, "type mismatch in array reference")

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
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! mode: f90
! license-default: "bsd"
! license-markup: "doxygen"
! End:
!
