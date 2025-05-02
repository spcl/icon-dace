# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/support/mo_util_uuid.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/build/verification//"
# 1 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/support/mo_util_uuid.f90"
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

! This module contains the routines for the calculation of 128bit
! data checksums/fingerprints (which are called "UUID's" in ICON).
!
! The larger part of this functionality is implemented in C routines
! in the "support" subdirectory, while this Fortran module acts
! merely as a wrapper. For *parallel* fingerprint calculation,
! however, the MPI-parallel communication is invoked on the Fortran
! level only.

MODULE mo_util_uuid

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_CHAR, C_SIGNED_CHAR, C_NULL_CHAR, &
    &                                    C_DOUBLE, C_INT, C_PTR, C_F_POINTER, C_LOC
  USE mo_util_uuid_types, ONLY: t_uuid
  USE mo_util_sort, ONLY: quicksort







  IMPLICIT NONE 

  PRIVATE

  PUBLIC :: uuid_generate
  PUBLIC :: compare_uuid
  PUBLIC :: uuid_parse
  PUBLIC :: uuid_unparse
  PUBLIC :: uuid2char
  PUBLIC :: char2uuid
  PUBLIC :: clear_uuid
  PUBLIC :: OPERATOR(==)
  PUBLIC :: UUID_EQUAL, UUID_EQUAL_LIMITED_ACCURACY, UUID_UNEQUAL

  !> module name string
  CHARACTER(LEN=*), PARAMETER :: modname = 'mo_util_uuid'

  TYPE, BIND(C) :: t_context
    INTEGER(C_INT) :: f64(2)
    INTEGER(C_INT) :: f48(2)
    INTEGER(C_INT) :: tl64(2)
    INTEGER(C_INT) :: tl48(2)
    INTEGER(C_INT) :: max_zero
  END TYPE t_context

  ENUM, BIND(c)
    ENUMERATOR :: UUID_EQUAL                  = 0
    ENUMERATOR :: UUID_EQUAL_LIMITED_ACCURACY = 1
    ENUMERATOR :: UUID_UNEQUAL                = 2
  END ENUM

  INTERFACE OPERATOR (==)
    MODULE PROCEDURE uuid_compare
  END INTERFACE OPERATOR (==)

  INTERFACE uuid_generate
    MODULE PROCEDURE uuid_generate_sequential
    MODULE PROCEDURE uuid_generate_parallel
    MODULE PROCEDURE uuid_generate_parallel1D
  END INTERFACE uuid_generate

  INTERFACE
    SUBROUTINE my_uuid_generate(val, nval, uuid) BIND(C,NAME='uuid_generate')
      IMPORT :: C_DOUBLE, C_INT, t_uuid
      REAL(c_double),   INTENT(IN)         :: val(*)
      INTEGER(c_int),   INTENT(IN), VALUE  :: nval
      TYPE(t_uuid),     INTENT(OUT)        :: uuid
    END SUBROUTINE my_uuid_generate
  END INTERFACE

  INTERFACE
    INTEGER(C_INT) FUNCTION my_compare_uuid(uuid_A, uuid_B, min_difference) BIND(C,NAME='compare_UUID')
      IMPORT :: t_uuid, C_INT, C_DOUBLE
      TYPE(t_uuid),     INTENT(IN), VALUE  :: uuid_A, uuid_B
      REAL(c_double),   INTENT(OUT)        :: min_difference
    END FUNCTION my_compare_uuid
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_unparse(uuid_string, uuid) BIND(C,NAME='uuid_unparse')
      IMPORT :: C_CHAR, t_uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(out) :: uuid_string
      TYPE(t_uuid),                    INTENT(in)  :: uuid
    END SUBROUTINE my_uuid_unparse
  END INTERFACE

  INTERFACE
    SUBROUTINE my_uuid_parse(uuid, uuid_string) BIND(C,NAME='uuid_parse')
      IMPORT :: C_CHAR, t_uuid
      TYPE(t_uuid),                    INTENT(out) :: uuid
      CHARACTER(C_CHAR), DIMENSION(*), INTENT(in)  :: uuid_string
    END SUBROUTINE my_uuid_parse
  END INTERFACE

  INTERFACE
    SUBROUTINE my_encode_uuid(fingerprint, uuid) BIND(C,NAME='encode_uuid')
      IMPORT :: t_uuid, C_PTR
      TYPE(C_PTR), VALUE        :: fingerprint
      TYPE(t_uuid), INTENT(out) :: uuid
    END SUBROUTINE my_encode_uuid
  END INTERFACE

  INTERFACE
    FUNCTION my_uuid_scan_data(val, nval) BIND(C,NAME='uuid_scan_data') RESULT(fingerprint)
      IMPORT :: C_PTR, C_DOUBLE, C_INT
      REAL(C_DOUBLE),     INTENT(IN) :: val(*)
      INTEGER(C_INT), INTENT(IN), VALUE :: nval
      TYPE(C_PTR) :: fingerprint
    END FUNCTION my_uuid_scan_data
  END INTERFACE

  INTERFACE
    FUNCTION my_concat_fingerprints(context0, context1) BIND(C,NAME='concat_fingerprints') RESULT(fingerprint)
      IMPORT :: C_PTR
      TYPE(C_PTR), VALUE :: context0
      TYPE(C_PTR), VALUE :: context1
      TYPE(C_PTR) :: fingerprint
    END FUNCTION my_concat_fingerprints
  END INTERFACE

  INTERFACE
    SUBROUTINE deallocate_c(ptr) BIND(C,NAME='deallocate_fingerprint') 
      IMPORT :: C_PTR
      TYPE(C_PTR), VALUE :: ptr
    END SUBROUTINE deallocate_c
  END INTERFACE

CONTAINS
 
  SUBROUTINE uuid_parse(uuid_string, uuid)
    CHARACTER(len=*), INTENT(in)  :: uuid_string
    TYPE(t_uuid),     INTENT(out) :: uuid
    CALL my_uuid_parse(uuid, TRIM(uuid_string)//C_NULL_CHAR)
  END SUBROUTINE uuid_parse

  SUBROUTINE uuid_unparse(uuid, uuid_string)
    TYPE(t_uuid),     INTENT(in)  :: uuid
    CHARACTER(len=*), INTENT(out) :: uuid_string
    CALL my_uuid_unparse(uuid_string, uuid)
  END SUBROUTINE uuid_unparse

  FUNCTION uuid_compare(uuid1, uuid2)
    LOGICAL :: uuid_compare
    TYPE(t_uuid),     INTENT(in)  :: uuid1
    TYPE(t_uuid),     INTENT(in)  :: uuid2
    uuid_compare = .TRUE.
!CDIR NOVECTOR
    IF (ANY(uuid1%data /= uuid2%data)) THEN
      uuid_compare = .FALSE.
    ENDIF
  END FUNCTION uuid_compare
  
  SUBROUTINE uuid2char(uuid, string)
    TYPE(t_uuid), INTENT(in) :: uuid
    CHARACTER(len=1), INTENT(out) :: string(16)
    string = TRANSFER(uuid%data, string)
  END SUBROUTINE uuid2char

  SUBROUTINE char2uuid(string, uuid)
    CHARACTER(len=1), INTENT(in) :: string(16)
    TYPE(t_uuid), INTENT(out) :: uuid
    uuid%data  = TRANSFER(string, uuid%data)
  END SUBROUTINE char2uuid

  SUBROUTINE uuid_generate_sequential(val, uuid)
    REAL(C_DOUBLE), INTENT(IN)  :: val(:)
    TYPE(t_uuid),   INTENT(out) :: uuid
    CALL my_uuid_generate(val, SIZE(val), uuid)
  END SUBROUTINE uuid_generate_sequential

  INTEGER(C_INT) FUNCTION compare_uuid(uuid_A, uuid_B) 
    TYPE(t_uuid),   INTENT(IN)  :: uuid_A, uuid_B    
    REAL(C_DOUBLE) :: min_difference
    compare_uuid = my_compare_uuid(uuid_A, uuid_B, min_difference) 
  END FUNCTION compare_uuid


  ! Parallel computation of UUID.
  !
  ! In contrast to the sequential version "uuid_generate", here the
  ! data may be spread over different PEs in an unordered way. Without
  ! using a global array, the UUID is then generated by concatenating
  ! independently computed fingerprints.
  !
  SUBROUTINE uuid_generate_parallel(comm, in_val, in_glbidx, glb_nval, uuid)
    INTEGER,        INTENT(IN)  :: comm             !< MPI communicator
    REAL(C_DOUBLE), INTENT(IN)  :: in_val(:,:)      !< input value, dim2: data-parallel
    INTEGER,        INTENT(IN)  :: in_glbidx(:)     !< "in_glbidx(:,i)": global index of "in_val(:,i)"
    INTEGER,        INTENT(IN)  :: glb_nval         !< total no. of global indices
    TYPE(t_uuid),   INTENT(OUT) :: uuid             !< (output:) UUID
    ! local variables
    INTEGER, PARAMETER :: ROOTPE    = 0
    INTEGER, PARAMETER :: dbg_level = 0

    CHARACTER(LEN=*), PARAMETER :: routine = modname//':uuid_generate_parallel'
    TYPE(C_PTR)                          :: fingerprint, fingerprint_merge, new_fingerprint
    TYPE(t_context), TARGET, ALLOCATABLE :: fingerprint_i(:)
    INTEGER                              :: nval, nrow, i, p_error, ipe, npes, target_pe,  &
      &                                     nval_local, chunksize, p_c_double_byte, p_c_double
    REAL(C_DOUBLE)                       :: tmp_c_double
    TYPE(t_context), POINTER             :: ptr
    INTEGER, ALLOCATABLE                 :: sendbuf(:), recvbuf(:,:), permutation(:),      &
      &                                     glbidx_sorted(:), glbidx_local(:), glbidx(:),  &
      &                                     sdispls(:), sendcounts(:), rdispls(:), recvcounts(:)
    REAL(C_DOUBLE), ALLOCATABLE          :: val_sorted(:,:), val_local(:,:)
    LOGICAL                              :: lcontiguous

    nrow = SIZE(in_val,1)
    nval = SIZE(in_val,2)

    IF (SIZE(in_glbidx) /= nval) THEN
      WRITE (0,*) routine, ": Size of input fields does not match!" ; RETURN
    END IF

# 443 "/home/primrose/Work/IconGrounds/icon-dace2/icon-model/support/mo_util_uuid.f90"
    ! non-MPI mode: execute sequential version
    ALLOCATE(permutation(nval))
    permutation(in_glbidx(:)) = (/ (i, i=1,nval) /)
    CALL uuid_generate_sequential(RESHAPE(in_val(:,permutation), (/ SIZE(in_val) /)), uuid)
    DEALLOCATE(permutation)

  END SUBROUTINE uuid_generate_parallel


  ! Parallel computation of UUID.
  ! Wrapper for 1D array input.
  !
  SUBROUTINE uuid_generate_parallel1D(comm, in_val, in_glbidx, glb_nval, uuid)
    INTEGER,        INTENT(IN)  :: comm
    REAL(C_DOUBLE), INTENT(IN)  :: in_val(:)
    INTEGER,        INTENT(IN)  :: in_glbidx(:)
    INTEGER,        INTENT(IN)  :: glb_nval
    TYPE(t_uuid),   INTENT(OUT) :: uuid
    ! local variables
    REAL(C_DOUBLE), ALLOCATABLE :: tmp_val(:,:)
    ALLOCATE(tmp_val(1,SIZE(in_val)))
    tmp_val(1,:) = in_val(:)
    CALL uuid_generate_parallel(comm, tmp_val, in_glbidx, glb_nval, uuid)
    DEALLOCATE(tmp_val)
  END SUBROUTINE uuid_generate_parallel1D


  SUBROUTINE clear_uuid(uuid)
    TYPE(t_uuid), INTENT(inout) :: uuid
    uuid%data(:) = INT(0, C_SIGNED_CHAR)
  END SUBROUTINE clear_uuid

END MODULE mo_util_uuid
