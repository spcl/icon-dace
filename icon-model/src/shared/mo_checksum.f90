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

! Utility functions to ensure that data does not change due to code modifications.

MODULE mo_checksum
    USE ISO_C_BINDING,  ONLY: C_INT32_T, C_DOUBLE, C_FLOAT, c_int, c_ptr, c_loc
    USE mo_mpi,         ONLY: p_comm_size, p_comm_rank, p_comm_work, p_gather
    USE mo_util_string, ONLY: int2string
    USE mo_cdi,         ONLY: DATATYPE_INT, DATATYPE_FLT64
#ifndef NOMPI
    USE mpi !, ONLY: MPI_INT32_T, MPI_GATHER
#endif

    IMPLICIT NONE

    PUBLIC printChecksum

    PRIVATE

    INTERFACE printChecksum
        MODULE PROCEDURE printChecksum_1d_int32
        MODULE PROCEDURE printChecksum_1d_float
        MODULE PROCEDURE printChecksum_1d_double

        MODULE PROCEDURE printChecksum_2d_int32
        MODULE PROCEDURE printChecksum_2d_float
        MODULE PROCEDURE printChecksum_2d_double

        MODULE PROCEDURE printChecksum_3d_int32
        MODULE PROCEDURE printChecksum_3d_float
        MODULE PROCEDURE printChecksum_3d_double

        MODULE PROCEDURE printChecksum_4d_int32
        MODULE PROCEDURE printChecksum_4d_float
        MODULE PROCEDURE printChecksum_4d_double

        MODULE PROCEDURE printChecksum_5d_int32
        MODULE PROCEDURE printChecksum_5d_float
        MODULE PROCEDURE printChecksum_5d_double
    END INTERFACE

    INTERFACE ! brain dead Fortran syntax for c bindings
      ! chksum is uint32_t in cdi, but this does not really matter for us.
      FUNCTION cdi_check_sum (cdi_type, cnt, buffer) BIND(c, name='cdiCheckSum') RESULT(chksum)
        IMPORT                :: c_ptr, c_int, c_int32_t
        INTEGER(c_int), VALUE :: cdi_type
        INTEGER(c_int), VALUE :: cnt
        TYPE(c_ptr),    VALUE :: buffer
        INTEGER(c_int32_t)    :: chksum
      END FUNCTION cdi_check_sum
    END INTERFACE

    CHARACTER(LEN = 1), PARAMETER :: kNibbles(0:15) = (/'0','1','2','3', &
      &                                                 '4','5','6','7', &
      &                                                 '8','9','a','b', &
      &                                                 'c','d','e','f'/)
    CHARACTER(LEN = *), PARAMETER :: moduleName = "mo_checksum"

CONTAINS

    SUBROUTINE checksum_to_string(chksum, str)
      INTEGER(c_int32_t), INTENT(in   ) :: chksum
      CHARACTER(len=8),   INTENT(  out) :: str
      INTEGER            :: i
      INTEGER(c_int32_t) :: chksum_

      chksum_ = chksum
      DO i = 1,8
        str(9-i:9-i) = kNibbles(IAND(15_c_int32_t, chksum_))
        chksum_      = ISHFT(chksum_, -4)
      END DO
    END SUBROUTINE

    ! If `opt_lDetails = .TRUE.` IS specified, this also prints a list of all the process local checksums.
    SUBROUTINE printChecksum_second_step(prefix, local_chksum, opt_comm, opt_lDetails)
        CHARACTER(LEN = *),        INTENT(in   ) :: prefix
        INTEGER(KIND = C_INT32_T), INTENT(in   ) :: local_chksum
        INTEGER,         OPTIONAL, INTENT(in   ) :: opt_comm
        LOGICAL,         OPTIONAL, INTENT(in   ) :: opt_lDetails

        CHARACTER(len=8)                               :: chksum_string
        INTEGER(KIND = C_INT32_T)                      :: hash
        INTEGER(KIND = C_INT32_T), ALLOCATABLE, TARGET :: processChecksums(:)
        INTEGER :: i, communicator, processCount
        LOGICAL :: printDetails

        !print a process local checksum
        printDetails = .FALSE.
        IF(PRESENT(opt_lDetails)) printDetails = opt_lDetails

#ifndef NOMPI
        !gather the process local checksums on process 0
        communicator = p_comm_work
        IF(PRESENT(opt_comm)) communicator = opt_comm
        processCount = p_comm_size(communicator)
        ALLOCATE(processChecksums(processCount))

        CALL p_gather(local_chksum, processChecksums, 0, communicator)

        !hash the results of the different processes down to a single VALUE AND print that.
        IF (p_comm_rank(communicator) == 0) THEN

          IF (ALL(processChecksums == 0)) THEN

            WRITE(0, *) prefix//' zero size array - no checksum available.'

          ELSE

            ! HACKHACKHACKHACK:: stupid serializeGetSizeInCore() from cdilib does not know about
            ! a DATATYPE_INT32, so use DATATYPE_INT as workaround
            ! (which SHOULD hopefully have same storage size...)
            hash = cdi_check_sum(DATATYPE_INT, processCount, c_loc(processChecksums))

            IF (printDetails) THEN
              WRITE(0, *) prefix//"details:"
              DO i = 1, processCount
                CALL checksum_to_string(processChecksums(i), chksum_string)
                WRITE(0, *) "checksum from process "//TRIM(int2string(i - 1))//": "//chksum_string
              END DO
            END IF

            !print the RESULT
            CALL checksum_to_string(hash, chksum_string)
            WRITE(0, *) prefix//chksum_string
          END IF
        ENDIF

        DEALLOCATE(processChecksums)
#else
        CALL checksum_to_string(local_chksum, chksum_string)
        WRITE(0, *) prefix//chksum_string
#endif
    END SUBROUTINE printChecksum_second_step


    ! AAAAAHHHHHH: Fortran! Many nice wrappers.
    SUBROUTINE printChecksum_int32(prefix, arr_size, array, opt_lDetails)
        CHARACTER(LEN = *),                INTENT(in   ) :: prefix
        INTEGER,                           INTENT(in   ) :: arr_size
        INTEGER(KIND = C_INT32_T), TARGET, INTENT(in   ) :: array(arr_size)
        LOGICAL,                 OPTIONAL, INTENT(in   ) :: opt_lDetails
        INTEGER(KIND = c_int32_t) :: local_chksum

        ! HACKHACKHACKHACK:: stupid serializeGetSizeInCore() from cdilib does
        ! not know about a DATATYPE_INT32, so use DATATYPE_INT as workaround
        ! (which SHOULD hopefully have same storage size...)
        IF (arr_size == 0) THEN
          local_chksum = 0
        ELSE
          local_chksum = cdi_check_sum(DATATYPE_INT, arr_size, c_loc(array))
        ENDIF
        CALL printChecksum_second_step(prefix, local_chksum, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_int32

    SUBROUTINE printChecksum_float(prefix, arr_size, array, opt_lDetails)
        CHARACTER(LEN = *),           INTENT(in   ) :: prefix
        INTEGER,                      INTENT(in   ) :: arr_size
        REAL(KIND = C_FLOAT), TARGET, INTENT(in   ) :: array(arr_size)
        LOGICAL,    OPTIONAL,         INTENT(in   ) :: opt_lDetails
        INTEGER(KIND = c_int32_t) :: local_chksum

        ! HACKHACKHACKHACK:: stupid serializeGetSizeInCore() from cdilib does
        ! not know about a DATATYPE_FLT32, so use DATATYPE_INT as workaround
        ! (which SHOULD hopefully have same storage size...)
        IF (arr_size == 0) THEN
          local_chksum = 0
        ELSE
          local_chksum = cdi_check_sum(DATATYPE_INT, arr_size, c_loc(array))
        ENDIF
        CALL printChecksum_second_step(prefix, local_chksum, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_float

    SUBROUTINE printChecksum_double(prefix, arr_size, array, opt_lDetails)
        CHARACTER(LEN = *),            INTENT(in   ) :: prefix
        INTEGER,                       INTENT(in   ) :: arr_size
        REAL(KIND = C_DOUBLE), TARGET, INTENT(in   ) :: array(arr_size)
        LOGICAL,     OPTIONAL,         INTENT(in   ) :: opt_lDetails
        INTEGER(KIND = c_int32_t) :: local_chksum

        IF (arr_size == 0) THEN
          local_chksum = 0
        ELSE
          local_chksum = cdi_check_sum(DATATYPE_FLT64, arr_size, c_loc(array))
        ENDIF
        CALL printChecksum_second_step(prefix, local_chksum, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_double

    ! This is the call layer computing the size from assumed shape arrays,
    ! which is then past down together with assumed size arrays.
    ! Fortran 2008s contiguous would simplify things, but coding standard...
    SUBROUTINE printChecksum_1d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),                INTENT(in   ) :: prefix
        INTEGER(KIND = C_INT32_T), TARGET, INTENT(in   ) :: array(:)
        LOGICAL,                 OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_int32(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_1d_int32

    SUBROUTINE printChecksum_1d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),   INTENT(in   ) :: prefix
        REAL(KIND = C_FLOAT), INTENT(in   ) :: array(:)
        LOGICAL,    OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_float(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_1d_float

    SUBROUTINE printChecksum_1d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),    INTENT(in   ) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(in   ) :: array(:)
        LOGICAL,     OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_double(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_1d_double


    SUBROUTINE printChecksum_2d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),                INTENT(in   ) :: prefix
        INTEGER(KIND = C_INT32_T), TARGET, INTENT(in   ) :: array(:,:)
        LOGICAL,                 OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_int32(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_2d_int32

    SUBROUTINE printChecksum_2d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),   INTENT(in   ) :: prefix
        REAL(KIND = C_FLOAT), INTENT(in   ) :: array(:,:)
        LOGICAL,    OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_float(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_2d_float

    SUBROUTINE printChecksum_2d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),    INTENT(in   ) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(in   ) :: array(:,:)
        LOGICAL,     OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_double(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_2d_double


    SUBROUTINE printChecksum_3d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),                INTENT(in   ) :: prefix
        INTEGER(KIND = C_INT32_T), TARGET, INTENT(in   ) :: array(:,:,:)
        LOGICAL,                 OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_int32(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_3d_int32

    SUBROUTINE printChecksum_3d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),   INTENT(in   ) :: prefix
        REAL(KIND = C_FLOAT), INTENT(in   ) :: array(:,:,:)
        LOGICAL,    OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_float(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_3d_float

    SUBROUTINE printChecksum_3d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),    INTENT(in   ) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(in   ) :: array(:,:,:)
        LOGICAL,     OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_double(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_3d_double


    SUBROUTINE printChecksum_4d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),                INTENT(in   ) :: prefix
        INTEGER(KIND = C_INT32_T), TARGET, INTENT(in   ) :: array(:,:,:,:)
        LOGICAL,                 OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_int32(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_4d_int32

    SUBROUTINE printChecksum_4d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),   INTENT(in   ) :: prefix
        REAL(KIND = C_FLOAT), INTENT(in   ) :: array(:,:,:,:)
        LOGICAL,    OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_float(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_4d_float

    SUBROUTINE printChecksum_4d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),    INTENT(in   ) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(in   ) :: array(:,:,:,:)
        LOGICAL,     OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_double(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_4d_double


    SUBROUTINE printChecksum_5d_int32(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),                INTENT(in   ) :: prefix
        INTEGER(KIND = C_INT32_T), TARGET, INTENT(in   ) :: array(:,:,:,:,:)
        LOGICAL,                 OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_int32(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_5d_int32

    SUBROUTINE printChecksum_5d_float(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),   INTENT(in   ) :: prefix
        REAL(KIND = C_FLOAT), INTENT(in   ) :: array(:,:,:,:,:)
        LOGICAL,    OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_float(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_5d_float

    SUBROUTINE printChecksum_5d_double(prefix, array, opt_lDetails)
        CHARACTER(LEN = *),    INTENT(in   ) :: prefix
        REAL(KIND = C_DOUBLE), INTENT(in   ) :: array(:,:,:,:,:)
        LOGICAL,     OPTIONAL, INTENT(in   ) :: opt_lDetails
        CALL printChecksum_double(prefix, size(array), array, opt_lDetails = opt_lDetails)
    END SUBROUTINE printChecksum_5d_double

END MODULE mo_checksum
