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

MODULE TEST_STRING
  USE FORTUTF
  USE mo_util_string

CONTAINS
  SUBROUTINE TEST_string_to_lower
    CHARACTER(len=10) :: uppercase, lowercase
    CALL TAG_TEST("TEST_to_lower")
    uppercase = 'ALLCAPITAL'
    lowercase = tolower(uppercase)
    CALL STRING_CONTAINS('allcapital', lowercase)
  END SUBROUTINE
  SUBROUTINE TEST_string_lowcase
    CHARACTER(len=10) :: testcase
    CALL TAG_TEST("TEST_lowcase")
    testcase = 'ALLCAPITAL'
    CALL lowcase(testcase)
    CALL STRING_CONTAINS('allcapital', testcase)
  END SUBROUTINE
  SUBROUTINE TEST_tocompact
    CHARACTER(len=30) :: testcase
    CALL TAG_TEST("TEST_tocompact")
    testcase = ' remove    extra   spaces '
    CALL tocompact(testcase)
    CALL STRING_CONTAINS('remove extra spaces', testcase)
  END SUBROUTINE
  SUBROUTINE TEST_int2string
    CHARACTER(len=6) :: testcase
    INTEGER :: n = 1234
    CALL TAG_TEST("TEST_int2string")
    testcase = int2string(n)
    CALL STRING_CONTAINS('1234', testcase)
  END SUBROUTINE
END MODULE TEST_STRING
