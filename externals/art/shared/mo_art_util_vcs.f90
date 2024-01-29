!
! mo_art_util_vcs
! This module provides the print out of vcs (e.g. git) information,
! like branch, hash and repository
! This module is based on src/shared/mo_util_vcs.f90 of ICON
! which is subject to the DWD and MPI-M-Software-License-Agreement in
! its most recent form.
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

MODULE mo_art_util_vcs

  USE, INTRINSIC :: iso_c_binding, ONLY: c_int, c_char, c_null_char

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: art_util_repository_url
  PUBLIC :: art_util_branch_name
  PUBLIC :: art_util_revision_key

  INTERFACE

  SUBROUTINE private_util_art_repository_url(name, actual_len) BIND(c,name='art_repository_url')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_art_repository_url

    SUBROUTINE private_util_art_branch_name(name, actual_len) BIND(c,name='art_branch_name')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_art_branch_name

    SUBROUTINE private_util_art_revision_key(name, actual_len) BIND(c,name='art_revision_key')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(inout) :: name
      INTEGER(c_int), INTENT(inout) :: actual_len
    END SUBROUTINE private_util_art_revision_key
  END INTERFACE

CONTAINS
    SUBROUTINE art_util_repository_url(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_art_repository_url(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE art_util_repository_url

  SUBROUTINE art_util_branch_name(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_art_branch_name(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
  END SUBROUTINE art_util_branch_name

  SUBROUTINE art_util_revision_key(name, actual_len)
    CHARACTER(len=*), INTENT(out) :: name
    INTEGER, INTENT(inout) :: actual_len
    INTEGER :: i
    CALL private_util_art_revision_key(name, actual_len)
    char_loop: DO i = 1 , LEN(name)
      IF (name(i:i) == c_null_char) EXIT char_loop
    END DO char_loop
    name(i:LEN(name)) = ' '
   END SUBROUTINE art_util_revision_key
END MODULE mo_art_util_vcs   
