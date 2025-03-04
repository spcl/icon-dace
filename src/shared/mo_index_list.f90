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






  IMPLICIT NONE

  PRIVATE

  PUBLIC :: generate_index_list, generate_index_list_batched

  INTERFACE generate_index_list




    MODULE PROCEDURE generate_index_list_i1_cpu
    MODULE PROCEDURE generate_index_list_i4_cpu

  END INTERFACE

! Warning: there will be no GPU -> CPU copy here, array of NUMBER
!  of indices will ONLY be on the GPU!!
  INTERFACE generate_index_list_batched




    MODULE PROCEDURE generate_index_list_batched_i1_cpu
    MODULE PROCEDURE generate_index_list_batched_i4_cpu

  END INTERFACE


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


END MODULE mo_index_list
