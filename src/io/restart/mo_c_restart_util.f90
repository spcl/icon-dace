! Interfaces to the functions defined in support/util_multifile_restart.c
!
! Initial implementation: Nathanael Huebbe
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

MODULE mo_c_restart_util
    USE ISO_C_BINDING, ONLY: C_CHAR, C_INT, C_NULL_CHAR

    IMPLICIT NONE

    PUBLIC :: createEmptyMultifileDir
    PUBLIC :: checkMultifileDir

    PRIVATE

    CHARACTER(*), PARAMETER :: modname = "mo_c_restart_util"

CONTAINS

    INTEGER FUNCTION createEmptyMultifileDir(path) RESULT(resultVar)
        CHARACTER(LEN = *), INTENT(IN) :: path

        CHARACTER(KIND = C_CHAR) :: pathCopy(LEN(path) + 1)
        INTEGER :: i

        INTERFACE
            INTEGER(C_INT) FUNCTION c_createEmptyMultifileDir(c_path) BIND(C, NAME = "createEmptyMultifileDir")
                IMPORT C_INT, C_CHAR
                CHARACTER(KIND = C_CHAR) :: c_path(*)
            END FUNCTION c_createEmptyMultifileDir
        END INTERFACE

        DO i = 1, LEN(path)
            pathCopy(i) = path(i:i)
        END DO
        pathCopy(LEN(path) + 1) = C_NULL_CHAR
        resultVar = c_createEmptyMultifileDir(pathCopy)
    END FUNCTION createEmptyMultifileDir

    INTEGER FUNCTION checkMultifileDir(path, expectedDomainCount, expectedFileCount) RESULT(resultVar)
        CHARACTER(LEN = *), INTENT(IN) :: path
        INTEGER, INTENT(IN) :: expectedDomainCount, expectedFileCount

        CHARACTER(KIND = C_CHAR) :: pathCopy(LEN(path) + 1)
        INTEGER :: i

        INTERFACE
            INTEGER(C_INT) FUNCTION c_checkMultifileDir(c_path, c_expectedDomainCount, c_expectedFileCount) &
              &    BIND(C, NAME = "checkMultifileDir")
                IMPORT C_INT, C_CHAR
                CHARACTER(KIND = C_CHAR) :: c_path(*)
                INTEGER, VALUE :: c_expectedDomainCount, c_expectedFileCount
            END FUNCTION c_checkMultifileDir
        END INTERFACE

        DO i = 1, LEN(path)
            pathCopy(i) = path(i:i)
        END DO
        pathCopy(LEN(path) + 1) = C_NULL_CHAR
        resultVar = c_checkMultifileDir(pathCopy, expectedDomainCount, expectedFileCount)
    END FUNCTION checkMultifileDir

END MODULE mo_c_restart_util
