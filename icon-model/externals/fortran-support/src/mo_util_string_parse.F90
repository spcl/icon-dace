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

!>
!! Utility module:
!!
!! Fortran-C-Interface for string parsing routine.
!!

MODULE mo_util_string_parse

  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int, c_char, c_null_char

  IMPLICIT NONE

  PRIVATE

  INTERFACE
    SUBROUTINE private_do_parse_intlist(parse_line, nvalues, nlev_val, out_values, ierr) &
      &  BIND(C, NAME='do_parse_intlist')
      IMPORT :: c_int, c_char
      CHARACTER(c_char), DIMENSION(*), INTENT(IN) :: parse_line
      INTEGER(c_int), VALUE, INTENT(IN)    :: nvalues
      INTEGER(c_int), VALUE, INTENT(IN)    :: nlev_val
      INTEGER(c_int), INTENT(INOUT) :: out_values(nvalues)
      INTEGER(c_int), INTENT(INOUT) :: ierr
    END SUBROUTINE private_do_parse_intlist
  END INTERFACE

  PUBLIC :: util_do_parse_intlist

CONTAINS

  ! ---------------------------------------------------------------------
  ! Subroutine parsing the string parse_line containing integer numbers.
  !
  !    Allowed is a comma- (or semicolon-) separated list of integers,
  !    and of integer ranges like "10...20".  One may also use the
  !    keyword "nlev" to denote the maximum integer (or, equivalently,
  !    "n" or "N").
  !
  !    Furthermore, arithmetic expressions like "(nlev - 2)" are
  !    possible.
  !
  !    Basic example:
  !       parse_line = "1,3,5...10,20...nlev"
  !    More complex example:
  !       parse_line = "1,2, 10 ...22;2;16-(3+11), N-2,16-(2+10);5"
  ! ---------------------------------------------------------------------
  SUBROUTINE util_do_parse_intlist(parse_line, nlev_value, out_values, ierr)
    CHARACTER(len=*), INTENT(IN)    :: parse_line
    INTEGER, INTENT(IN)    :: nlev_value !< number to substitute for "N"/"nlev"
    INTEGER, INTENT(INOUT) :: out_values(0:)
    INTEGER, INTENT(INOUT) :: ierr
    CALL private_do_parse_intlist(TRIM(parse_line)//c_null_char, SIZE(out_values), &
      &                           nlev_value, out_values, ierr)
  END SUBROUTINE util_do_parse_intlist

END MODULE mo_util_string_parse
