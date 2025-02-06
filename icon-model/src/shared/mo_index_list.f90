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

MODULE mo_index_list

  USE mo_kind, ONLY: i1, i2, i4

#ifdef _OPENACC
  USE openacc
  USE, INTRINSIC :: iso_c_binding
#endif

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: generate_index_list, generate_index_list_batched

  INTERFACE generate_index_list
#ifdef _OPENACC
    MODULE PROCEDURE generate_index_list_i1
    MODULE PROCEDURE generate_index_list_i4
#else
    MODULE PROCEDURE generate_index_list_i1_cpu
    MODULE PROCEDURE generate_index_list_i4_cpu
#endif
  END INTERFACE

! Warning: there will be no GPU -> CPU copy here, array of NUMBER
!  of indices will ONLY be on the GPU!!
  INTERFACE generate_index_list_batched
#ifdef _OPENACC
    MODULE PROCEDURE generate_index_list_batched_i1
    MODULE PROCEDURE generate_index_list_batched_i4
#else
    MODULE PROCEDURE generate_index_list_batched_i1_cpu
    MODULE PROCEDURE generate_index_list_batched_i4_cpu
#endif
  END INTERFACE

#ifdef _OPENACC

  INTERFACE generate_index_list_gpu

    SUBROUTINE generate_index_list_gpu_single                  &
        & (conditions, startid, endid,                          &
        &  indices, nvalid, data_size,                          &
        &  copy_to_host, stream)                                &
        & BIND(C, name="c_generate_index_list_gpu_single")

      USE openacc
      USE iso_c_binding

      TYPE(c_ptr),              INTENT(in),  VALUE  :: conditions
      INTEGER(c_int),           INTENT(in),  VALUE  :: startid
      INTEGER(c_int),           INTENT(in),  VALUE  :: endid
      TYPE(c_ptr),              INTENT(in),  VALUE  :: indices ! actually intent is more like inout, but this doesn't work now
      TYPE(c_ptr),              INTENT(in),  VALUE  :: nvalid
      INTEGER(c_int),           INTENT(in),  VALUE  :: data_size
      LOGICAL(c_bool),          INTENT(in),  VALUE  :: copy_to_host
      INTEGER(acc_handle_kind), INTENT(in),  VALUE  :: stream

    END SUBROUTINE generate_index_list_gpu_single

    SUBROUTINE generate_index_list_gpu_batched                  &
        & (batch_size,                                           &
        &  conditions, cond_stride,                              &
        &  startid, endid,                                       &
        &  indices, idx_stride,                                  &
        &  nvalid, data_size,                                    & 
        &  stream )                                              &
        & BIND(C, name="c_generate_index_list_gpu_batched")

      USE openacc
      USE iso_c_binding

      INTEGER(c_int),           INTENT(in),  VALUE  :: batch_size
      TYPE(c_ptr),              INTENT(in),  VALUE  :: conditions
      INTEGER(c_int),           INTENT(in),  VALUE  :: cond_stride
      INTEGER(c_int),           INTENT(in),  VALUE  :: startid
      INTEGER(c_int),           INTENT(in),  VALUE  :: endid
      TYPE(c_ptr),              INTENT(in),  VALUE  :: indices ! actually intent is more like inout, but this doesn't work now
      INTEGER(c_int),           INTENT(in),  VALUE  :: idx_stride
      TYPE(c_ptr),              INTENT(in),  VALUE  :: nvalid
      INTEGER(c_int),           INTENT(in),  VALUE  :: data_size
      INTEGER(acc_handle_kind), INTENT(in),  VALUE  :: stream

    END SUBROUTINE generate_index_list_gpu_batched
  END INTERFACE

#endif

  CONTAINS


! Regular CPU implementation with a simple loop

  SUBROUTINE generate_index_list_i1_cpu(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_acc_copy_to_host, opt_use_acc)
    INTEGER(i1), INTENT(in)           :: conditions(:)
    INTEGER,     INTENT(inout)        :: indices(:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid
    ! These arguments are used in the OpenACC variant, but not in the CPU one 
    INTEGER,     INTENT(in), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(in), OPTIONAL :: opt_acc_copy_to_host
    LOGICAL,     INTENT(in), OPTIONAL :: opt_use_acc

    INTEGER :: i, nvalid_loc

    nvalid_loc = 0
    DO i = startid, endid
      IF (conditions(i) /= 0) THEN
        nvalid_loc = nvalid_loc + 1
        indices(nvalid_loc) = i
      END IF
    END DO
    nvalid = nvalid_loc
  END SUBROUTINE generate_index_list_i1_cpu

  SUBROUTINE generate_index_list_i4_cpu(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_acc_copy_to_host, opt_use_acc)
    INTEGER(i4), INTENT(in)           :: conditions(:)
    INTEGER,     INTENT(inout)        :: indices(:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid
    ! These arguments are used in the OpenACC variant, but not in the CPU one 
    INTEGER,     INTENT(in), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(in), OPTIONAL :: opt_acc_copy_to_host
    LOGICAL,     INTENT(in), OPTIONAL :: opt_use_acc

    INTEGER :: i, nvalid_loc

    nvalid_loc = 0
    DO i = startid, endid
      IF (conditions(i) /= 0) THEN
        nvalid_loc = nvalid_loc + 1
        indices(nvalid_loc) = i
      END IF
    END DO
    nvalid = nvalid_loc
  END SUBROUTINE generate_index_list_i4_cpu

  SUBROUTINE generate_index_list_batched_i1_cpu(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_use_acc)
    INTEGER(i1), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(inout)        :: nvalid(:)
    ! This argument is used in the OpenACC variant, but not in the GPU one 
    INTEGER,     INTENT(in), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(IN), OPTIONAL :: opt_use_acc

    INTEGER :: i, batch, batch_size
    batch_size = size(conditions, 2)
    nvalid = 0

    DO batch = 1, batch_size
      CALL generate_index_list_i1_cpu(              &
        conditions(:,batch), indices(:, batch), &
        startid, endid, nvalid(batch) )
    END DO

  END SUBROUTINE generate_index_list_batched_i1_cpu

  SUBROUTINE generate_index_list_batched_i4_cpu(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_use_acc)
    INTEGER(i4), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(inout)        :: nvalid(:)
    ! This argument is used in the OpenACC variant, but not in the GPU one
    INTEGER,     INTENT(in), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(IN), OPTIONAL :: opt_use_acc

    INTEGER :: i, batch, batch_size
    batch_size = size(conditions, 2)
    nvalid = 0

    DO batch = 1, batch_size
      CALL generate_index_list_i4_cpu(              &
        conditions(:,batch), indices(:, batch), &
        startid, endid, nvalid(batch) )
    END DO

  END SUBROUTINE generate_index_list_batched_i4_cpu

#ifdef _OPENACC

  ! on the gpu call the cub library through c++

  ! 1 byte
  SUBROUTINE generate_index_list_i1(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_acc_copy_to_host, opt_use_acc)
    INTEGER(i1), INTENT(IN)           :: conditions(:)
    INTEGER,     INTENT(INOUT)        :: indices(:)
    INTEGER,     INTENT(IN)           :: startid
    INTEGER,     INTENT(IN)           :: endid
    INTEGER,     INTENT(OUT)          :: nvalid
    INTEGER,     INTENT(IN), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(IN), OPTIONAL :: opt_use_acc
    LOGICAL,     INTENT(in), OPTIONAL :: opt_acc_copy_to_host

    INTEGER(acc_handle_kind) :: stream
    LOGICAL :: use_acc, acc_copy_to_host

    IF (PRESENT(opt_use_acc)) THEN
      use_acc = opt_use_acc
    ELSE
      use_acc = .TRUE.
    ENDIF

    ! run on GPU
    IF (use_acc) THEN

#if !defined( __HIP__ )
      IF ( PRESENT(opt_acc_async_queue) ) THEN
        stream = acc_get_cuda_stream(opt_acc_async_queue)
      ELSE
        stream = acc_get_cuda_stream(acc_async_sync)
      END IF
#endif

      IF ( PRESENT(opt_acc_copy_to_host) ) THEN
        acc_copy_to_host = opt_acc_copy_to_host
      ELSE
        acc_copy_to_host = .TRUE.
      END IF

      !$ACC HOST_DATA USE_DEVICE(conditions, indices)
      !$ACC HOST_DATA USE_DEVICE(nvalid) IF(.not. acc_copy_to_host)
      CALL generate_index_list_gpu_single(       &
        & c_loc(conditions(1)),                   &
        & startid, endid,                         &
        & c_loc(indices(1)), c_loc(nvalid),       &
        & 1, logical(acc_copy_to_host, c_bool),   &
        & stream)
      !$ACC END HOST_DATA
      !$ACC END HOST_DATA

#ifdef __HIP__ 
      CALL acc_wait_all()
#endif
    ELSE ! run on CPU
      CALL generate_index_list_i1_cpu(conditions, indices, startid, endid, nvalid)
    ENDIF

  END SUBROUTINE generate_index_list_i1

  ! 4 bytes
  SUBROUTINE generate_index_list_i4(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_acc_copy_to_host, opt_use_acc)
    INTEGER(i4), INTENT(in)           :: conditions(:)
    INTEGER,     INTENT(inout)        :: indices(:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid
    INTEGER,     INTENT(in), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(in), OPTIONAL :: opt_acc_copy_to_host
    LOGICAL,     INTENT(in), OPTIONAL :: opt_use_acc

    INTEGER(acc_handle_kind) :: stream
    LOGICAL :: use_acc, acc_copy_to_host

    IF (PRESENT(opt_use_acc)) THEN
      use_acc = opt_use_acc
    ELSE
      use_acc = .TRUE.
    ENDIF

    ! run on GPU
    IF (use_acc) THEN

#if !defined( __HIP__ )
      IF ( PRESENT(opt_acc_async_queue) ) THEN
        stream = acc_get_cuda_stream(opt_acc_async_queue)
      ELSE
        stream = acc_get_cuda_stream(acc_async_sync)
      END IF
#endif

      IF ( PRESENT(opt_acc_copy_to_host) ) THEN
        acc_copy_to_host = opt_acc_copy_to_host
      ELSE
        acc_copy_to_host = .TRUE.
      END IF

      !$ACC HOST_DATA USE_DEVICE(conditions, indices)
      !$ACC HOST_DATA USE_DEVICE(nvalid) IF(.not. acc_copy_to_host)
      CALL generate_index_list_gpu_single(       &
        & c_loc(conditions),                      &
        & startid, endid,                         &
        & c_loc(indices), c_loc(nvalid),          &
        & 4, logical(acc_copy_to_host, c_bool),   &
        & stream)
      !$ACC END HOST_DATA
      !$ACC END HOST_DATA

#ifdef __HIP__ 
      CALL acc_wait_all()
#endif
    ELSE ! run on CPU
      CALL generate_index_list_i4_cpu(conditions, indices, startid, endid, nvalid)
    ENDIF
  END SUBROUTINE generate_index_list_i4

  !
  ! Batched
  !

  ! 1 byte
  SUBROUTINE generate_index_list_batched_i1(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_use_acc)
    ! Attention: all the values, including nvalid (!!)
    !  are updated asynchronously on device ONLY
    INTEGER(i1), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid(:)
    INTEGER,     INTENT(in), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(in), OPTIONAL :: opt_use_acc

    INTEGER(acc_handle_kind) :: stream
    INTEGER :: i
    INTEGER :: batch_size
    INTEGER :: cond_stride, idx_stride
    LOGICAL :: use_acc

    IF (PRESENT(opt_use_acc)) THEN
      use_acc = opt_use_acc
    ELSE
      use_acc = .TRUE.
    ENDIF

    batch_size = size(conditions, 2)

    ! run on GPU
    IF (use_acc) THEN

#if !defined( __HIP__ )
      IF ( PRESENT(opt_acc_async_queue) ) THEN
        stream = acc_get_cuda_stream(opt_acc_async_queue)
      ELSE
        stream = acc_get_cuda_stream(acc_async_sync)
      END IF
#endif

      ! Hacky way to support non-contiguous slices
      IF ( batch_size > 1 ) THEN
        cond_stride = ( LOC(conditions(1,2)) - LOC(conditions(1,1)) ) / SIZEOF(conditions(1,1))
        idx_stride  = ( LOC(indices   (1,2)) - LOC(indices   (1,1)) ) / SIZEOF(indices   (1,1))
      END IF

      !$ACC HOST_DATA USE_DEVICE(conditions, indices, nvalid)
      CALL generate_index_list_gpu_batched(         &
          & batch_size,                              &
          & c_loc(conditions(1,1)), cond_stride,     &
          & startid, endid,                          &
          & c_loc(indices(1,1)), idx_stride,         &
          & c_loc(nvalid(1)), 1, stream )
      !$ACC END HOST_DATA

#ifdef __HIP__ 
      CALL acc_wait_all()
#endif
    ELSE ! run on CPU
      CALL generate_index_list_batched_i1_cpu(conditions, indices, startid, endid, nvalid)
    ENDIF
  END SUBROUTINE generate_index_list_batched_i1

  ! 4 bytes
  SUBROUTINE generate_index_list_batched_i4(conditions, indices, startid, endid, nvalid, opt_acc_async_queue, opt_use_acc)
    ! Attention: all the values, including nvalid (!!)
    !  are updated asynchronously on device ONLY
    INTEGER(i4), INTENT(in)           :: conditions(:,:)
    INTEGER,     INTENT(inout)        :: indices(:,:)
    INTEGER,     INTENT(in)           :: startid
    INTEGER,     INTENT(in)           :: endid
    INTEGER,     INTENT(out)          :: nvalid(:)
    INTEGER,     INTENT(in), OPTIONAL :: opt_acc_async_queue
    LOGICAL,     INTENT(in), OPTIONAL :: opt_use_acc

    INTEGER(acc_handle_kind) :: stream
    INTEGER :: i
    INTEGER :: batch_size
    INTEGER :: cond_stride, idx_stride
    LOGICAL :: use_acc

    IF (PRESENT(opt_use_acc)) THEN
      use_acc = opt_use_acc
    ELSE
      use_acc = .TRUE.
    ENDIF

    batch_size = size(conditions, 2)

    ! run on GPU
    IF (use_acc) THEN

#if !defined( __HIP__ )
      IF ( PRESENT(opt_acc_async_queue) ) THEN
        stream = acc_get_cuda_stream(opt_acc_async_queue)
      ELSE
        stream = acc_get_cuda_stream(acc_async_sync)
      END IF
#endif

      ! Hacky way to support non-contiguous slices
      IF ( batch_size > 1 ) THEN
        cond_stride = ( LOC(conditions(1,2)) - LOC(conditions(1,1)) ) / SIZEOF(conditions(1,1))
        idx_stride  = ( LOC(indices   (1,2)) - LOC(indices   (1,1)) ) / SIZEOF(indices   (1,1))
      END IF

      !$ACC HOST_DATA USE_DEVICE(conditions, indices, nvalid)
      CALL generate_index_list_gpu_batched(         &
          & batch_size,                              &
          & c_loc(conditions(1,1)), cond_stride,     &
          & startid, endid,                          &
          & c_loc(indices(1,1)), idx_stride,         &
          & c_loc(nvalid(1)), 4, stream )
      !$ACC END HOST_DATA

#ifdef __HIP__ 
      CALL acc_wait_all()
#endif
    ELSE ! run on CPU
      CALL generate_index_list_batched_i4_cpu(conditions, indices, startid, endid, nvalid)
    ENDIF
  END SUBROUTINE generate_index_list_batched_i4

#endif

END MODULE mo_index_list
