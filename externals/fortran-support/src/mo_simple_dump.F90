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

MODULE mo_simple_dump

  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: wp => real64

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: dump2text

  INTERFACE dump2text
    MODULE PROCEDURE dump2text1_wp
    MODULE PROCEDURE dump2text2_wp
    MODULE PROCEDURE dump2text3_wp
  END INTERFACE dump2text

CONTAINS

  SUBROUTINE dump2text1_wp(array, name, preserve_original)
    REAL(wp), INTENT(inout) :: array(:)
    CHARACTER(len=*), INTENT(in) :: name
    LOGICAL, INTENT(in), OPTIONAL :: preserve_original

    LOGICAL :: do_copy
    REAL(wp), ALLOCATABLE :: backup(:)

    do_copy = .TRUE.
    IF (PRESENT(preserve_original)) do_copy = preserve_original

    IF (do_copy) backup = array

    !$ACC UPDATE HOST(array)
    OPEN (17, file=name)
    WRITE (17, *) array
    CLOSE (17)

    IF (do_copy) array = backup
  END SUBROUTINE dump2text1_wp

  SUBROUTINE dump2text2_wp(array, name, preserve_original)
    REAL(wp), INTENT(inout) :: array(:, :)
    CHARACTER(len=*), INTENT(in) :: name
    LOGICAL, INTENT(in), OPTIONAL :: preserve_original

    LOGICAL :: do_copy
    REAL(wp), ALLOCATABLE :: backup(:, :)

    do_copy = .TRUE.
    IF (PRESENT(preserve_original)) do_copy = preserve_original

    IF (do_copy) backup = array

    !$ACC UPDATE HOST(array)
    OPEN (17, file=name)
    WRITE (17, *) array
    CLOSE (17)

    IF (do_copy) array = backup
  END SUBROUTINE dump2text2_wp

  SUBROUTINE dump2text3_wp(array, name, preserve_original)
    REAL(wp), INTENT(inout) :: array(:, :, :)
    CHARACTER(len=*), INTENT(in) :: name
    LOGICAL, INTENT(in), OPTIONAL :: preserve_original

    LOGICAL :: do_copy
    REAL(wp), ALLOCATABLE :: backup(:, :, :)

    do_copy = .TRUE.
    IF (PRESENT(preserve_original)) do_copy = preserve_original

    IF (do_copy) backup = array

    !$ACC UPDATE HOST(array)
    OPEN (17, file=name)
    WRITE (17, *) array
    CLOSE (17)

    IF (do_copy) array = backup
  END SUBROUTINE dump2text3_wp
END MODULE mo_simple_dump
